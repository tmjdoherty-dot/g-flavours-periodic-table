"""Average-of-configuration Hartree-Fock for atoms with g electron spin flavours.

Scaled units (s = g r, E/g^2, z = Z/g): one-electron part has nuclear charge z, every
electron-electron term carries gamma = 1/g. Shell a = (n, l) holds q_a electrons (q_a <= g(2l+1)).

E = sum_a q_a I_a + gamma [ sum_a q_a(q_a-1)/2 E_aa + sum_{a<b} q_a q_b E_ab ]
E_aa = A_l [ g F0(aa) - sum_k (l k l;000)^2 F^k(aa) ],  A_l = (2l+1)/(g(2l+1)-1)   (0 if capacity 1)
E_ab = F0(ab) - (1/g) sum_k (la k lb;000)^2 G^k(ab)
(g = 2 gives Slater's average-of-configuration energy.) Fractional q allowed (same polynomial).

Orbitals live on the log grid of atom.py in the w-representation (P^2 dr per point = w^2).
Minimisation: Newton-like preconditioned steps in the eigenbasis of the LDA Hamiltonian
(occupied-virtual), explicit rotations for occupied-occupied pairs with different q,
Lowdin re-orthonormalisation, Armijo backtracking on E. Optional fractional q: projected
gradient on q with dE/dq.
"""
import numpy as np
from math import factorial as fac
from scipy.linalg import eigh_tridiagonal, eigh, solve_banded
from atom import Grid, scf as lda_scf, solve_levels, project

LN = "spdfg"


def tj0(l1, k, l2):
    """(l1 k l2; 0 0 0)^2"""
    J = l1 + k + l2
    if J % 2 or k < abs(l1 - l2) or k > l1 + l2:
        return 0.0
    g = J // 2
    return fac(J - 2 * l1) * fac(J - 2 * k) * fac(J - 2 * l2) / fac(J + 1) * (
        fac(g) / (fac(g - l1) * fac(g - k) * fac(g - l2))) ** 2


def kin_tridiag(grid, l):
    r, h = grid.r, grid.h
    d = (0.5 * (2 / h**2 + (l + 0.5) ** 2)) / r**2
    e = (-0.5 / h**2) / (r[:-1] * r[1:])
    return d, e


def tri_mv(d, e, X):
    """tridiagonal (d, e) times X (N x m)"""
    Y = d[:, None] * X
    Y[:-1] += e[:, None] * X[1:]
    Y[1:] += e[:, None] * X[:-1]
    return Y


class Yk:
    def __init__(self, r):
        self.r = r
        self.pw = {}

    def __call__(self, rho, k):
        """rho: (P, N) pair densities (per point). Returns Y^k (P, N)."""
        if k not in self.pw:
            self.pw[k] = (self.r ** k, self.r ** (-k - 1.0))
        rk, rmk = self.pw[k]
        inner = np.cumsum(rho * rk, axis=1) * rmk
        t = rho * rmk
        outer = np.cumsum(t[:, ::-1], axis=1)[:, ::-1] - t
        return inner + outer * rk


class Atom:
    def __init__(self, z, g, shells, grid=None, lda=None, qf=None):
        """shells: list of (n, l, q). Shells with q = 0 are kept as (virtual) orbitals.
        qf: optional (S, g) flavour occupations -> single flavour-resolved configuration
        (Hund mode, e.g. maximal flavour polarisation) instead of the configuration average."""
        self.qf = None if qf is None else np.asarray(qf, float)
        self.z, self.g, self.gam = z, g, 1.0 / g
        self.grid = grid or Grid(N=1200)
        r = self.grid.r
        self.r = r
        self.shells = [(n, l) for n, l, q in shells]
        self.q = np.array([q for n, l, q in shells], float)
        self.cap = np.array([g * (2 * l + 1) for n, l in self.shells], float)
        self.ls = sorted(set(l for n, l in self.shells))
        self.idx = {l: [i for i, (n, ll) in enumerate(self.shells) if ll == l] for l in self.ls}
        for l in self.ls:  # order by n inside each l
            self.idx[l].sort(key=lambda i: self.shells[i][0])
        self.Y = Yk(r)
        self.A = np.array([(2 * l + 1) / (g * (2 * l + 1) - 1) if g * (2 * l + 1) > 1 else 0.0
                           for n, l in self.shells])
        self.Bl = np.array([(2 * l + 1) / (2 * l) if l > 0 else 0.0 for n, l in self.shells])
        # reference LDA (lam = 1/g) for guess + preconditioner
        if lda == "fixed":  # LDA with this configuration's occupations (per flavour) held fixed
            lv = solve_levels(self.grid, -z / r)
            f0 = np.zeros(len(lv))
            pos = {(n_, l_): i for i, (n_, l_, e_, w_) in enumerate(lv)}
            for (n_, l_), qq in zip(self.shells, self.q):
                f0[pos[(n_, l_)]] = qq / g
            lda = lda_scf(z, 1.0 / g, self.grid, f0=f0, fixed=True, mix=0.2, tol=1e-6, maxit=800)
        elif lda is None:
            lda = lda_scf(z, 1.0 / g, self.grid, eta=8, tol=1e-6, maxit=600)
        V = lda["V"]
        self.kin = {l: kin_tridiag(self.grid, l) for l in self.ls}
        self.U, self.e = {}, {}
        W = np.zeros((len(r), len(self.shells)))
        for l in self.ls:
            d, e = self.kin[l]
            self.U[l], self.e[l] = (d + V, e), None
            nmax = max(self.shells[i][0] for i in self.idx[l])
            k = nmax - l
            vals, vecs = eigh_tridiagonal(d + V, e, select="i", select_range=(0, k - 1),
                                          lapack_driver="stebz", tol=1e-13)
            for i in self.idx[l]:
                W[:, i] = vecs[:, self.shells[i][0] - l - 1]
        self.W = W
        # k lists
        self.kself = [[k for k in range(0, 2 * l + 1, 2)] for n, l in self.shells]
        S = len(self.shells)
        self.pairs = []
        for a in range(S):
            for b in range(a + 1, S):
                la, lb = self.shells[a][1], self.shells[b][1]
                ks = [k for k in range(abs(la - lb), la + lb + 1, 2)]
                self.pairs.append((a, b, ks, [tj0(la, k, lb) for k in ks]))
        self.cself = [[tj0(l, k, l) for k in ks] for (n, l), ks in zip(self.shells, self.kself)]

    # ------------------------------------------------------------------
    def evaluate(self, W, q):
        g, gam = self.g, self.gam
        S = W.shape[1]
        rho_d = (W * W).T  # (S, N) shell densities
        # self Y's
        Yself = {}
        kmax = max(max(ks) for ks in self.kself)
        for k in range(0, kmax + 1, 2):
            sel = [a for a in range(S) if k in self.kself[a]]
            Yk = self.Y(rho_d[sel], k)
            for j, a in enumerate(sel):
                Yself[(a, k)] = Yk[j]
        Y0 = np.array([Yself[(a, 0)] for a in range(S)])  # (S, N)
        # pair exchange Y's grouped by k
        bucket = {}
        for p, (a, b, ks, ds) in enumerate(self.pairs):
            for k in ks:
                bucket.setdefault(k, []).append(p)
        Yab = {}
        for k, plist in bucket.items():
            rho = np.array([W[:, self.pairs[p][0]] * W[:, self.pairs[p][1]] for p in plist])
            Yk = self.Y(rho, k)
            for j, p in enumerate(plist):
                Yab[(p, k)] = Yk[j]
        # one-electron
        HW = np.empty_like(W)
        for l in self.ls:
            d, e = self.kin[l]
            ii = self.idx[l]
            HW[:, ii] = tri_mv(d - self.z / self.r, e, W[:, ii])
        I = np.einsum("ia,ia->a", W, HW)
        # two-electron energies
        F0 = rho_d @ Y0.T  # F0[a,b] = int rho_a Y0_b
        Xs = [sum(c * Yself[(a, k)] for k, c in zip(self.kself[a], self.cself[a])) for a in range(S)]
        Fk = np.array([rho_d[a] @ Xs[a] for a in range(S)])  # sum_k c_k F^k(aa)
        Gab = np.zeros((S, S))
        Xab = {}
        for p, (a, b, ks, ds) in enumerate(self.pairs):
            X = sum(dk * Yab[(p, k)] for k, dk in zip(ks, ds))
            Xab[p] = X
            Gab[a, b] = Gab[b, a] = (W[:, a] * W[:, b]) @ X
        # coefficients: E_intra = al*F0aa - be*Fk ; E_inter = q_a q_b F0ab - Smat_ab Gab
        if self.qf is None:
            pairs_ = q * (q - 1) / 2
            al, be = pairs_ * self.A * g, pairs_ * self.A
            al_q, be_q = (q - 1) * self.A * g / 2, (q - 1) * self.A / 2  # al/q, be/q (limits ok at q=0)
            Smat = np.outer(q, q) / g
            Sq = np.tile(q, (S, 1)) / g  # Sq[a,b] = Smat[a,b]/q_a
        else:
            qf = self.qf
            Ps = np.sum(qf * (qf - 1) / 2, axis=1)
            Pd = (q ** 2 - np.sum(qf ** 2, axis=1)) / 2
            al, be = Ps * self.Bl + Pd, Ps * self.Bl
            qs = np.where(q > 0, q, 1.0)
            al_q, be_q = al / qs, be / qs
            Smat = qf @ qf.T
            Sq = Smat / qs[:, None]
        Fdiag = np.diag(F0).copy()
        F0o = F0.copy()
        np.fill_diagonal(F0o, 0.0)
        np.fill_diagonal(Gab, 0.0)
        E = q @ I + gam * (np.sum(al * Fdiag - be * Fk) + 0.5 * (q @ F0o @ q) - 0.5 * np.sum(Smat * Gab))
        Eaa = self.A * (g * Fdiag - Fk)
        Eab = F0o - Gab / g
        dEdq = I + gam * ((q - 0.5) * Eaa + Eab @ q)  # valid for the configuration average
        VH = q @ Y0
        FW = HW.copy()
        for a in range(S):
            selfpot = 2 * al_q[a] * Y0[a] - 2 * be_q[a] * Xs[a]
            FW[:, a] += gam * (selfpot + VH - q[a] * Y0[a]) * W[:, a]
        for p, (a, b, ks, ds) in enumerate(self.pairs):
            FW[:, a] -= gam * Sq[a, b] * Xab[p] * W[:, b]
            FW[:, b] -= gam * Sq[b, a] * Xab[p] * W[:, a]
        return E, FW, dEdq

    # ------------------------------------------------------------------
    def step_dirs(self, W, FW, q):
        """Preconditioned descent directions. Returns D (ov part), Theta per l, grad norm."""
        D = np.zeros_like(W)
        thetas = {}
        gmax = 0.0
        for l in self.ls:
            ii = self.idx[l]
            C = W[:, ii]
            R = FW[:, ii]
            M = C.T @ R  # M[b,a] = <b|F_a|a>
            eps = np.diag(M).copy()
            Rp = R - C @ M
            dd, ee = self.U[l]
            Dl = np.empty_like(Rp)
            for j in range(len(ii)):
                # Olsen / Jacobi-Davidson correction: x = -(H-e)^-1 R + (H-e)^-1 C alpha, x orth C
                ab = np.zeros((3, len(dd)))
                ab[0, 1:] = ee
                ab[1] = dd - eps[j]
                ab[2, :-1] = ee
                Y = solve_banded((1, 1), ab, np.column_stack([Rp[:, j], C]))
                y, Z = Y[:, 0], Y[:, 1:]
                al = np.linalg.solve(C.T @ Z, C.T @ y)
                Dl[:, j] = -(y - Z @ al)
                if np.dot(Dl[:, j], Rp[:, j]) >= 0:  # indefinite reference: fall back to gradient
                    Dl[:, j] = -Rp[:, j]
            dec = np.abs(np.einsum("ia,ia->a", Rp, Dl))
            gmax = max(gmax, np.sqrt(dec.max()))
            D[:, ii] = Dl
            n = len(ii)
            Th = np.zeros((n, n))
            for x in range(n):
                for y in range(x + 1, n):
                    qa, qb = q[ii[x]], q[ii[y]]
                    if abs(qa - qb) < 1e-10 or (qa == 0 and qb == 0):
                        continue
                    # rotation a' = a cos t + b sin t: dE/dt = 2(qa<b|Fa|a> - qb<a|Fb|b>)
                    gab = 2 * (qa * M[y, x] - qb * M[x, y])
                    hab = 2 * abs((qa - qb) * (eps[y] - eps[x])) + 1e-3
                    t = np.clip(-gab / hab, -0.2, 0.2)
                    gmax = max(gmax, abs(gab) / max(qa, qb, 1e-12))
                    Th[y, x] = t
                    Th[x, y] = -t
            thetas[l] = Th
        return D, thetas, gmax

    def apply(self, W, D, thetas, alpha):
        Wn = W.copy()
        for l in self.ls:
            ii = self.idx[l]
            C = W[:, ii] + alpha * D[:, ii]
            Th = alpha * thetas[l]
            if np.any(Th):
                ev, V = eigh(1j * Th)  # expm of skew matrix
                R = (V @ np.diag(np.exp(-1j * ev)) @ V.conj().T).real
                C = C @ R
            S = C.T @ C
            sv, sV = eigh(S)
            C = C @ (sV @ np.diag(sv ** -0.5) @ sV.T)
            Wn[:, ii] = C
        return Wn

    def solve(self, free=None, tol=1e-6, maxit=3000, eta=1.0, verbose=False):
        """Minimise E. free: indices of shells whose q is optimised (fractional), or None."""
        W, q = self.W.copy(), self.q.copy()
        E, FW, dEdq = self.evaluate(W, q)
        hist = []
        a_last = 0.5
        for it in range(maxit):
            D, th, gmax = self.step_dirs(W, FW, q)
            qg = 0.0
            if free is not None:
                fq = q[free]
                part = (fq > 1e-9) & (fq < self.cap[free] - 1e-9)
                dq = dEdq[free]
                if part.any():
                    qg = dq[part].max() - dq[part].min()
            if gmax < tol and qg < tol:
                break
            if (len(hist) > 30 and abs(hist[-30] - E) < 1e-12 * max(1.0, abs(E)) and qg < tol
                    and gmax < 1e-4):
                break
            alpha, ok = min(1.0, 2 * a_last), False
            for _ in range(14):
                Wn = self.apply(W, D, th, alpha)
                En, FWn, dEdqn = self.evaluate(Wn, q)
                if En <= E + 1e-13:
                    ok = True
                    break
                alpha *= 0.5
            if ok:
                W, E, FW, dEdq = Wn, En, FWn, dEdqn
                a_last = alpha
            else:
                a_last = 1.0
            if free is not None:
                cap = self.cap[free]
                tot = q[free].sum()
                st = eta
                for _ in range(12):
                    qn = q.copy()
                    qn[free] = project(q[free] - st * dEdq[free], cap, tot)
                    En, FWn, dEdqn = self.evaluate(W, qn)
                    if En <= E + 1e-13:
                        q, E, FW, dEdq = qn, En, FWn, dEdqn
                        eta = min(st * 1.5, 50.0)
                        break
                    st *= 0.5
                else:
                    eta = max(st, 1e-3)
            if verbose and it % 50 == 0:
                print(it, E, gmax, qg, alpha)
            hist.append(E)
        self.W, self.q = W, q
        eps = np.einsum("ia,ia->a", W, FW)
        return dict(E=E, q=q, eps=eps, dEdq=dEdq, it=it, gmax=gmax, conv=bool(gmax < 1e-4))


def hund_qf(shells, g):
    """maximal flavour polarisation: each shell fills flavour 1 to 2l+1, then flavour 2, ..."""
    qf = np.zeros((len(shells), g))
    for i, (n, l, q) in enumerate(shells):
        left = q
        for s in range(g):
            take = min(2 * l + 1, left)
            qf[i, s] = take
            left -= take
    return qf


def core(ncore_per_flavour_shells, g):
    """helper: full shells list -> (n,l,q) with q = capacity"""
    return [(n, l, g * (2 * l + 1)) for n, l in ncore_per_flavour_shells]


AR = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1)]
KR = AR + [(3, 2), (4, 0), (4, 1)]
XE = KR + [(4, 2), (5, 0), (5, 1)]
