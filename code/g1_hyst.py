"""G1b: hysteresis at the 3d collapse. Sweep z up and down through the 3d onset with continuation."""
import time, numpy as np
from concurrent.futures import ProcessPoolExecutor
from atom import Grid, scf
LN = "spdfg"

def sweep(args):
    lam, zs = args
    g = Grid(N=1200); V = f = None; rows = []
    for z in zs:
        o = scf(z, lam, g, V0=V, f0=f, eta=4, mix=0.25, tol=1e-6, maxit=1500)
        V, f = o["V"], o["f"]
        rows.append((z, o["occ"][(3, 2)][0], o["occ"][(4, 1)][0], o["occ"][(4, 0)][0], o["occ"][(3, 2)][1], o["err"]))
    return lam, zs[0] < zs[-1], rows

if __name__ == "__main__":
    t = time.time()
    up = np.round(np.arange(9.75, 11.76, 0.0833333), 4)
    jobs = [(lam, s) for lam in (0.0, 0.25, 0.5) for s in (up, up[::-1])]
    with ProcessPoolExecutor(6) as ex:
        res = list(ex.map(sweep, jobs))
    lines = []
    for lam in (0.0, 0.25, 0.5):
        U = {r[0]: r for l, u, rows in res if l == lam and u for r in rows}
        Dn = {r[0]: r for l, u, rows in res if l == lam and not u for r in rows}
        lines.append(f"lam={lam}:  z   3d_up 4p_up | 3d_down 4p_down   (eps3d up/down, err)")
        for z in sorted(U):
            a, b = U[z], Dn[z]
            flag = "  <-- differs" if abs(a[1] - b[1]) > 0.05 else ""
            lines.append(f"  {z:6.3f}  {a[1]:.2f} {a[2]:.2f} | {b[1]:.2f} {b[2]:.2f}   ({a[4]:+.4f}/{b[4]:+.4f}, {max(a[5], b[5]):.0e}){flag}")
    txt = "\n".join(lines); print(txt)
    open("results/g1_hysteresis.txt", "w").write(txt)
    print("time", time.time() - t)
