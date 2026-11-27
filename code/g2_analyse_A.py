"""G2 part A: ground configurations (n-1)d^x ns^y per g from HF (Hund and configuration average).
Energies extrapolated from N = 1600, 2400 (error ~ h^2). Exception = ground s occupancy differs from
Madelung (s fills first: y = min(delta, g))."""
import json, sys
from collections import defaultdict
import numpy as np

FILES = sys.argv[1:] or ["results/g2/A.json"]
R, ext, conv_ok, src = [], defaultdict(dict), defaultdict(lambda: True), defaultdict(dict)
nlower = 0
for fn in FILES:
    Rf = json.load(open(fn))
    R += Rf
    E = defaultdict(dict)
    for r in Rf:
        if "err" in r:
            continue
        E[(r["g"], r["row"], r["delta"], r["mode"], r["x"], r["y"])][r["N"]] = (r["E"], r["conv"])
    for (g, row, d, mode, x, y), v in E.items():
        if 1600 in v and 2400 in v and v[1600][1] and v[2400][1]:
            e1, e2 = v[1600][0], v[2400][0]
            ex = e2 + (e2 - e1) / 1.25
            k = (g, row, d, mode)
            if y not in ext[k] or ex < ext[k][y] - 1e-6:
                if y in ext[k]:
                    nlower += 1
                ext[k][y] = ex
                src[k][y] = fn
bad = [r for r in R if "err" in r or not r.get("conv", False)]
NAMES = {"3d": "K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn".split(), "4d": "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd".split()}
EXP_Y = {"3d": [1, 2, 2, 2, 2, 1, 2, 2, 2, 2, 1, 2], "4d": [1, 2, 2, 2, 1, 1, 2, 1, 1, 0, 1, 2]}

out = []
P = lambda s="": (print(s), out.append(s))
P(f"files: {FILES}; rows: {len(R)}, failed/unconverged: {len(bad)}; later file lower by >1e-6: {nlower}")
gs = sorted({k[0] for k in ext})
for mode in ("hund", "avg"):
    P(f"\n=== {mode} ===")
    for row in ("3d", "4d"):
        P(f"-- {row} row: s-occupancy of the ground configuration per element (Madelung = s full first)")
        for g in gs:
            ys, exc, gaps = [], [], []
            for d in range(1, 6 * g + 1):
                k = (g, row, d, mode)
                if k not in ext or not ext[k]:
                    ys.append("?")
                    continue
                yy = min(ext[k], key=ext[k].get)
                ym = min(d, g)
                ys.append(str(yy) if yy < 10 else "*")
                if yy != ym and ym in ext[k]:
                    exc.append(d)
                    gaps.append((ext[k][ym] - ext[k][yy]) * g * g * 27.2114)  # eV, real units
            nd = 5 * g  # d-block members (delta = g+1 .. 6g)
            nexc_d = sum(1 for d in exc if d > g)
            P(f"  g={g}: y = {' '.join(ys)}")
            P(f"        exceptions at delta = {exc}  ({nexc_d}/{nd} of the d block = {nexc_d/nd:.2f})"
              + (f", max gain {max(gaps):.2f} eV" if gaps else ""))
            if g == 2 and all(ext.get((2, row, d, mode)) for d in range(1, 13)):
                yg = [min(ext[(2, row, d, mode)], key=ext[(2, row, d, mode)].get) for d in range(1, 13)]
                hit = [NAMES[row][i] for i in range(12) if yg[i] == EXP_Y[row][i]]
                miss = [f"{NAMES[row][i]}(HF {yg[i]}, exp {EXP_Y[row][i]})" for i in range(12) if yg[i] != EXP_Y[row][i]]
                P(f"        vs experiment: {len(hit)}/12 s-occupancies right; wrong: {', '.join(miss) or 'none'}")
out_path = "results/g2/A_summary.txt"
open(out_path, "w", encoding="utf-8").write("\n".join(out))
