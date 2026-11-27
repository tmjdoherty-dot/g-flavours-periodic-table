"""G2 part B: fractional (configuration-average) HF scans. Per g and row:
 - outer-s drain: max over the d block of (g - q_s)/g  (1 = s empties)   [compare LDA result 2]
 - p head start: max q_p/g before the d shell holds 1 electron per flavour  [compare LDA result 4]
 - d onset: first z with q_d/g > 0.1 that stays > 0.1 for the next 0.5   [compare LDA result 1]"""
import json, sys
from collections import defaultdict
import numpy as np

R = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "results/g2/B.json"))
err = [r for r in R if "err" in r]
D = defaultdict(list)
for r in R:
    if "err" not in r:
        D[(r["g"], r["row"])].append(r)
core = {"3d": 9, "4d": 18}
out = []
P = lambda s="": (print(s), out.append(s))
P(f"points: {len(R)}, errors: {len(err)}, unconverged: {sum(1 for r in R if 'err' not in r and not r['conv'])}")
for row in ("3d", "4d"):
    P(f"\n-- {row} row (per flavour; LDA g=2 drain was 0.58 for 3d, 1.00 for 4d)")
    P("   g   s-drain  p-headstart  d-onset(z-core)  flagged")
    for g in sorted({k[0] for k in D}):
        rs = sorted(D[(g, row)], key=lambda r: r["z"])
        z = np.array([r["z"] for r in rs]) - core[row]
        qs = np.array([r["qs"] for r in rs]) / g
        qd = np.array([r["qd"] for r in rs]) / g
        qp = np.array([r["qp"] for r in rs]) / g
        ok = np.array([r["conv"] for r in rs])
        # robust: 5-point running median over all points (isolated unconverged spikes removed)
        qsm = np.array([np.median(qs[max(0, i - 2):i + 3]) for i in range(len(qs))])
        blk = (z > 1.2) & (z < 5.8)
        drain = (1 - qsm[blk]).max() if blk.any() else np.nan
        qpm = np.array([np.median(qp[max(0, i - 2):i + 3]) for i in range(len(qp))])
        pre = (qd < 1.0)
        head = qpm[pre].max() if pre.any() else np.nan
        on = np.nan
        for i in range(len(z)):
            w = (z >= z[i]) & (z <= z[i] + 0.5)
            if qd[i] > 0.1 and (qd[w] > 0.1).all():
                on = z[i]
                break
        P(f"  {g:2d}   {drain:6.2f}   {head:8.2f}     {on:7.3f}       {(~ok).sum()}")
p = (sys.argv[1] if len(sys.argv) > 1 else "results/g2/B.json").replace(".json", "_summary.txt")
open(p, "w", encoding="utf-8").write("\n".join(out))
