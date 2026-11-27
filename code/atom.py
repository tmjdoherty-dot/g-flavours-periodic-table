"""Spherical mean-field atom with g electron spin flavours.

Per-flavour units (s = g r, E/g^2, z = Z/g): the g-flavour atom with Slater (LDA) exchange
is the g = 1 problem with nuclear charge z and exchange multiplied by lam = 1/g:

    [-1/2 d^2 + l(l+1)/2s^2 - z/s + V_H[rho] - lam (6/pi)^(1/3) rho^(1/3)] u = eps u

rho = per-flavour density (integrates to z), shell (n,l) holds up to 2l+1 per flavour.
lam = 0 is the g -> infinity limit (pure Hartree, exact mean field), lam = 1/2 is real
chemistry (exchange only, no correlation unless corr=True, which is only for the g = 2 check).

Radial equation on a log grid x = ln s, u = s^(1/2) v, w = s v:
    tridiagonal symmetric matrix, eigenvalues eps directly (Dirichlet at both ends).
Occupations: Fermi smearing (kT small) over all (n,l) levels; fractional at the Fermi level.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal

LMAX = 4
NMAX = 9


class Grid:
    def __init__(self, rmin=1e-6, rmax=120.0, N=1800):
        self.x = np.linspace(np.log(rmin), np.log(rmax), N)
        self.h = self.x[1] - self.x[0]
        self.r = np.exp(self.x)
        self.N = N


def pz81_corr(n):
    """Perdew-Zunger 1981 unpolarized correlation potential (Hartree). Only for g = 2 check."""
    n = np.maximum(n, 1e-30)
    rs = (3.0 / (4 * np.pi * n)) ** (1 / 3)
    g, b1, b2 = -0.1423, 1.0529, 0.3334
    A, B, C, D = 0.0311, -0.048, 0.0020, -0.0116
    v = np.empty_like(rs)
    hi = rs >= 1
    sq = np.sqrt(rs[hi])
    den = 1 + b1 * sq + b2 * rs[hi]
    ec = g / den
    v[hi] = ec * (1 + 7 / 6 * b1 * sq + 4 / 3 * b2 * rs[hi]) / den
    lo = ~hi
    lr = np.log(rs[lo])
    v[lo] = A * lr + (B - A / 3) + 2 / 3 * C * rs[lo] * lr + (2 * D - C) / 3 * rs[lo]
    return v


def hartree(grid, dq):
    """dq = 4 pi r^2 rho dr per grid point (sums to electron count). Returns V_H(r)."""
    q = np.cumsum(dq)  # enclosed charge (left Riemann; fine at this resolution)
    outer = np.cumsum((dq / grid.r)[::-1])[::-1]
    return q / grid.r + outer - dq / grid.r


def solve_levels(grid, V):
    """Lowest levels for l = 0..LMAX. Returns list of (n, l, eps, w)."""
    r, h = grid.r, grid.h
    out = []
    off_base = -0.5 / h**2
    for l in range(LMAX + 1):
        d = (0.5 * (2 / h**2 + (l + 0.5) ** 2)) / r**2 + V + 0.0
        e = off_base / (r[:-1] * r[1:])
        k = NMAX - l
        vals, vecs = eigh_tridiagonal(d, e, select="i", select_range=(0, k - 1),
                                       lapack_driver="stebz", tol=1e-13)
        for i in range(k):
            out.append((l + 1 + i, l, vals[i], vecs[:, i]))
    return out


def occupy(levels, nel, kT):
    eps = np.array([lv[2] for lv in levels])
    cap = np.array([2 * lv[1] + 1 for lv in levels], float)
    lo, hi = eps.min() - 50 * kT - 1, eps.max() + 50 * kT + 1
    for _ in range(200):
        mu = 0.5 * (lo + hi)
        f = cap / (1 + np.exp(np.clip((eps - mu) / kT, -700, 700)))
        if f.sum() > nel:
            hi = mu
        else:
            lo = mu
    return f, mu


def project(y, cap, nel):
    """Euclidean projection onto {sum f = nel, 0 <= f <= cap}."""
    lo, hi = (y - cap).min() - 1, y.max() + 1
    for _ in range(100):
        t = 0.5 * (lo + hi)
        s = np.clip(y - t, 0, cap).sum()
        if s > nel:
            lo = t
        else:
            hi = t
    return np.clip(y - 0.5 * (lo + hi), 0, cap)


def scf(z, lam, grid, V0=None, f0=None, kT=1e-3, g_corr=None, tol=1e-7, maxit=600,
        mix=0.35, eta=0.5, fixed=False):
    """Self-consistent neutral atom, minimising E over fractional occupations f
    (Janak: dE/df_i = eps_i) by projected gradient steps, joint with density mixing.
    Result = fractional-aufbau ground state: shells below mu full, shells at mu share it.
    g_corr: int -> add PZ81 correlation (only for the g = 2 NIST check)."""
    r = grid.r
    Vn = -z / r
    if V0 is None:
        V = -z / r * np.exp(-r * max(z, 0.5) ** (1 / 3))
        V = np.minimum(V, -1e-12)
    else:
        V = V0.copy()
    f = None if f0 is None else f0.copy()
    dq_old = None
    err = np.nan
    for it in range(maxit):
        levels = solve_levels(grid, V)
        eps = np.array([lv[2] for lv in levels])
        cap = np.array([2 * lv[1] + 1 for lv in levels], float)
        if f is None:
            f, _ = occupy(levels, z, kT)
            df = 1.0
        elif fixed:
            df = 0.0
        else:
            fn = project(f - eta * eps, cap, z)
            df = np.abs(fn - f).max()
            f = fn
        dq = np.zeros_like(r)
        for fi, lv in zip(f, levels):
            if fi > 1e-14:
                dq += fi * lv[3] ** 2
        if dq_old is not None:
            err = np.abs(dq - dq_old).sum()
            dq = mix * dq + (1 - mix) * dq_old
        rho = dq / (4 * np.pi * r**3 * grid.h)
        Vnew = Vn + hartree(grid, dq) - lam * (6 / np.pi) ** (1 / 3) * np.cbrt(rho)
        if g_corr is not None:
            g = g_corr  # real n = g^4 rho_scaled; real potential / g^2 -> scaled
            Vnew = Vnew + pz81_corr(g**4 * rho) / g**2
        V = Vnew
        dq_old = dq
        if err < tol and df < 1e-8:
            break
    part = (f > 1e-6) & (f < cap - 1e-6)
    mu = eps[part].mean() if part.any() else eps[f > 1e-6].max()
    occ = {(lv[0], lv[1]): (fi, lv[2]) for fi, lv in zip(f, levels)}
    return dict(V=V, f=f, occ=occ, mu=mu, it=it, err=err, df=df)


LNAME = "spdfg"


def config_str(occ, thr=0.02):
    items = sorted(((n, l, f) for (n, l), (f, e) in occ.items() if f > thr), key=lambda t: (t[0], t[1]))
    return " ".join(f"{n}{LNAME[l]}{f:.2f}" for n, l, f in items)
