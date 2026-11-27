"""G2 on Modal: Hartree-Fock version of the g-flavour periodic table.

A: integer configurations [core] (n-1)d^x ns^y, all y = 0..g, Hund (max flavour polarisation) and
   configuration-average energies, two grids (N = 1600, 2400) for Richardson extrapolation.
   Rows: 3d (Ar core) and 4d (Kr core); g in {1,2,3,4,5,6,8}.
B: fractional (configuration-average) HF scans through the 3d and 4d rows: free shells ns, (n-1)d, np.
   g in {1,2,3,4,6,8,16}; z step 1/24; N = 1600.
C: lanthanide row, Hund: [Xe] 6s^g 4f^x 5d^y, y in {0,1,2}; g in {1,2,3,4,6}; two grids.
Run from chemistry/:  python -m modal run modal_jobs/hf_modal.py
"""
import json, os, sys, time
import modal

app = modal.App("idea-g2-hf")
image = (modal.Image.debian_slim(python_version="3.11")
         .pip_install("numpy==1.26.4", "scipy==1.13.1")
         .env({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
         .add_local_python_source("atom", "hf"))

AR = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1)]
KR = AR + [(3, 2), (4, 0), (4, 1)]
XE = KR + [(4, 2), (5, 0), (5, 1)]
ROWS = {"3d": (AR, (4, 0), (3, 2), (4, 1)), "4d": (KR, (5, 0), (4, 2), (5, 1))}


def cpf(shells):  # per-flavour electrons in full shells
    return sum(2 * l + 1 for n, l in shells)


@app.function(image=image, cpu=1.0, memory=1024, timeout=3600, retries=1)
def task_A(g, row, delta, guess=None):
    """one element: delta = electrons beyond the core (Z = g*cpf(core) + delta)."""
    import numpy as np
    from atom import Grid
    from hf import Atom, core, hund_qf
    cr, s_sh, d_sh, p_sh = ROWS[row]
    Z = g * cpf(cr) + delta
    z = Z / g
    out = []
    for y in range(0, g + 1):
        x = delta - y
        if x < 0 or x > 5 * g:
            continue
        sh = core(cr, g) + [(d_sh[0], d_sh[1], x), (s_sh[0], s_sh[1], y)]
        sh = [t for t in sh if t[2] > 0]
        for mode in ("hund", "avg"):
            for N in (1600, 2400):
                t = time.time()
                qf = hund_qf(sh, g) if mode == "hund" else None
                try:
                    a = Atom(z, g, sh, grid=Grid(rmin=1e-9, N=N), qf=qf, lda=guess)
                    o = a.solve(tol=1e-7, maxit=1500)
                    out.append(dict(g=g, row=row, delta=delta, Z=Z, x=x, y=y, mode=mode, N=N,
                                    E=float(o["E"]), conv=bool(o["conv"]), gmax=float(o["gmax"]),
                                    it=int(o["it"]), t=time.time() - t))
                except Exception as ex:
                    out.append(dict(g=g, row=row, delta=delta, Z=Z, x=x, y=y, mode=mode, N=N,
                                    err=repr(ex)))
    return out


@app.function(image=image, cpu=1.0, memory=1024, timeout=3600, retries=1)
def task_B(g, row, z):
    import numpy as np
    from atom import Grid
    from hf import Atom, core
    cr, s_sh, d_sh, p_sh = ROWS[row]
    delta = g * (z - cpf(cr))  # electrons beyond the core (may be fractional)
    qs = min(delta, g)
    qd = min(max(delta - g, 0.0), 5 * g)
    qp = max(delta - 6 * g, 0.0)
    sh = core(cr, g) + [(s_sh[0], s_sh[1], qs), (d_sh[0], d_sh[1], qd), (p_sh[0], p_sh[1], qp)]
    nc = len(cr)
    t = time.time()
    try:
        a = Atom(z, g, sh, grid=Grid(rmin=1e-9, N=1600))
        o = a.solve(free=[nc, nc + 1, nc + 2], tol=1e-7, maxit=4000)
        return dict(g=g, row=row, z=z, qs=float(o["q"][nc]), qd=float(o["q"][nc + 1]),
                    qp=float(o["q"][nc + 2]), eps=[float(v) for v in o["eps"][nc:nc + 3]],
                    dEdq=[float(v) for v in o["dEdq"][nc:nc + 3]], E=float(o["E"]),
                    conv=bool(o["conv"]), it=int(o["it"]), t=time.time() - t)
    except Exception as ex:
        return dict(g=g, row=row, z=z, err=repr(ex))


@app.function(image=image, cpu=1.0, memory=1024, timeout=3600, retries=1)
def task_C(g, delta):
    """[Xe] 6s^g 4f^x 5d^y, delta = x + y."""
    from atom import Grid
    from hf import Atom, core, hund_qf
    Z = g * (cpf(XE) + 1) + delta
    z = Z / g
    out = []
    for y in (0, 1, 2):
        x = delta - y
        if x < 0 or x > 7 * g or y > 5 * g:
            continue
        sh = core(XE, g) + [(6, 0, g), (4, 3, x), (5, 2, y)]
        sh = [t for t in sh if t[2] > 0]
        for N in (1600, 2400):
            t = time.time()
            try:
                a = Atom(z, g, sh, grid=Grid(rmin=1e-9, N=N), qf=hund_qf(sh, g), lda="fixed")
                o = a.solve(tol=1e-7, maxit=1500)
                out.append(dict(g=g, delta=delta, Z=Z, x=x, y=y, N=N, E=float(o["E"]),
                                conv=bool(o["conv"]), gmax=float(o["gmax"]), it=int(o["it"]),
                                t=time.time() - t))
            except Exception as ex:
                out.append(dict(g=g, delta=delta, Z=Z, x=x, y=y, N=N, err=repr(ex)))
    return out


@app.local_entrypoint()
def main(parts: str = "ABC", test: bool = False, guess: str = ""):
    t0 = time.time()
    os.makedirs("results/g2", exist_ok=True)
    if "A" in parts:
        gs = [2] if test else [1, 2, 3, 4, 5, 6, 8]
        argsA = [(g, row, d, guess or None) for g in gs for row in (["3d"] if test else ["3d", "4d"])
                 for d in range(1, 6 * g + 1)]
        if test:
            argsA = [a for a in argsA if a[2] in (6, 11)]
    # simple blocking runs (starmap streams results)
    res = {}
    if "A" in parts:
        res["A"] = [r for rows in task_A.starmap(argsA, order_outputs=False) for r in rows]
        json.dump(res["A"], open("results/g2/A%s%s.json" % ("_test" if test else "", "_" + guess if guess else ""), "w"))
        print(f"A done: {len(res['A'])} rows, {time.time()-t0:.0f}s", flush=True)
    if "B" in parts:
        import numpy as np
        gs = [2] if test else [1, 2, 3, 4, 6, 8, 16]
        zs3 = [round(9 + k / 24, 6) for k in range(0, 157)]
        zs4 = [round(18 + k / 24, 6) for k in range(0, 157)]
        argsB = [(g, "3d", z) for g in gs for z in zs3] + [(g, "4d", z) for g in gs for z in zs4]
        if test:
            argsB = argsB[::40]
        res["B"] = list(task_B.starmap(argsB, order_outputs=False))
        json.dump(res["B"], open("results/g2/B%s.json" % ("_test" if test else ""), "w"))
        print(f"B done: {len(res['B'])} rows, {time.time()-t0:.0f}s", flush=True)
    if "C" in parts:
        gs = [2] if test else [1, 2, 3, 4, 6]
        argsC = [(g, d) for g in gs for d in range(1, 7 * g + 1)]
        if test:
            argsC = argsC[:3]
        res["C"] = [r for rows in task_C.starmap(argsC, order_outputs=False) for r in rows]
        json.dump(res["C"], open("results/g2/C%s.json" % ("_test" if test else ""), "w"))
        print(f"C done: {len(res['C'])} rows, {time.time()-t0:.0f}s", flush=True)
    print("total", time.time() - t0)
