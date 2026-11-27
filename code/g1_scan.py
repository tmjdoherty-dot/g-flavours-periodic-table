"""G1: periodic table with g spin flavours. Scan z = Z/g for lam = 1/g in {0, 1/4, 1/3, 1/2, 1},
plus the real g = 2 table with PZ81 correlation (integer Z). 12 worker processes, each runs one
(lam, z-segment) with continuation in z. Output: results/g1_scan.npz"""
import time, sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from atom import Grid, scf, solve_levels

LAMS = {"inf": 0.0, "4": 0.25, "3": 1 / 3, "2": 0.5, "1": 1.0}
ZS = np.round(np.concatenate([np.arange(3, 361) / 12, np.arange(180, 361) / 6]), 6)  # 1/12 to z=30, 1/6 above
SEG = 15


def run_one(z, lam, grid, V, f, corr):
    o = scf(z, lam, grid, V0=V, f0=None, eta=8, tol=1e-5, maxit=400, g_corr=corr)
    if not (o["err"] < 1e-5):  # gentler retry; still-unconverged points are kept and flagged
        o = scf(z, lam, grid, V0=o["V"], f0=o["f"], eta=4, mix=0.15, tol=1e-5, maxit=400, g_corr=corr)
    return o


def task(args):
    key, lam, zs, corr, fn = args
    grid = Grid(N=1200)
    V = f = None
    rows = []
    for z in zs:
        t = time.time()
        o = run_one(z, lam, grid, V, f, corr)
        V = o["V"]
        f = o["f"]
        eps = np.array([v[1] for v in o["occ"].values()])
        rows.append((key, z, o["f"].copy(), eps, o["mu"], o["err"], o["it"], time.time() - t))
    np.savez(fn, key=np.array([r[0] for r in rows]), z=np.array([r[1] for r in rows]),
             f=np.array([r[2] for r in rows]), eps=np.array([r[3] for r in rows]),
             mu=np.array([r[4] for r in rows]), err=np.array([r[5] for r in rows]),
             it=np.array([r[6] for r in rows]), t=np.array([r[7] for r in rows]))
    return fn, rows[0][1], rows[-1][1], sum(r[7] for r in rows)


def labels():
    g = Grid(N=200)
    return [(lv[0], lv[1]) for lv in solve_levels(g, -1 / g.r)]


if __name__ == "__main__":
    import os, glob
    from concurrent.futures import as_completed
    os.makedirs("results/g1_parts", exist_ok=True)
    t0 = time.time()
    jobs = []
    for key, lam in LAMS.items():
        for i in range(0, len(ZS), SEG):
            jobs.append((key, lam, ZS[i:i + SEG], None, f"results/g1_parts/{key}_{i:04d}.npz"))
    zc = np.arange(1, 121) / 2.0
    for i in range(0, len(zc), 10):
        jobs.append(("2c", 0.5, zc[i:i + 10], 2, f"results/g1_parts/2c_{i:04d}.npz"))
    jobs = [j for j in jobs if not os.path.exists(j[4])]
    # slow ones first: lam = 0 at high z, then by z
    jobs.sort(key=lambda j: -(j[2].mean() * (4 if j[1] == 0 else 1)))
    print(f"{len(jobs)} jobs", flush=True)
    with ProcessPoolExecutor(max_workers=12) as ex:
        futs = [ex.submit(task, j) for j in jobs]
        for n, fu in enumerate(as_completed(futs)):
            fn, za, zb, cpu = fu.result()
            print(f"[{n+1}/{len(jobs)}] {fn} z {za:.2f}-{zb:.2f} cpu {cpu:.0f}s  wall {time.time()-t0:.0f}s", flush=True)
    parts = [np.load(p) for p in sorted(glob.glob("results/g1_parts/*.npz"))]
    cat = {k: np.concatenate([p[k] for p in parts]) for k in ("key", "z", "f", "eps", "mu", "err", "it", "t")}
    np.savez("results/g1_scan.npz", labels=np.array(labels()), **cat)
    print("total", time.time() - t0)
