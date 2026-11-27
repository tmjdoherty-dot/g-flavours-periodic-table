"""G1 analysis v2: persistent onsets (f > 0.1 and stays > 0.1 for the next 0.5 in z),
s-shell drain across each d/f row (how many s electrons per flavour leave), converged points only."""
import numpy as np
import sys
FN = sys.argv[1] if len(sys.argv) > 1 else "results/g1_scan_run2.npz"
d = np.load(FN, allow_pickle=True)
LN = "spdfg"
lab = [f"{n}{LN[l]}" for n, l in d["labels"]]
keys = ["inf", "4", "3", "2", "2c", "1"]
D = {}
for k in keys:
    m = d["key"] == k; o = np.argsort(d["z"][m])
    D[k] = {q: d[q][m][o] for q in ("z", "f", "err", "mu")}
out = []
P = lambda s="": (print(s), out.append(s))

def onset(k, s, thr=0.1, hold=0.5):
    z, f = D[k]["z"], D[k]["f"][:, lab.index(s)]
    for i in range(len(z)):
        w = (z >= z[i]) & (z <= z[i] + hold)
        if f[i] > thr and (f[w] > thr).all():
            return z[i]
    return np.nan

MAD = {"2p": 2, "3d": 10, "4p": 15, "4d": 19, "5p": 24, "4f": 28, "5d": 35, "6p": 40, "5f": 44, "6d": 51, "7p": 56, "5g": 60}
P("Persistent onsets (per-flavour z; real element Z = g*z):")
P("shell  Mad " + " ".join(f"{k:>6}" for k in keys))
for s, mz in MAD.items():
    P(f"{s:>5} {mz:4d} " + " ".join(f"{onset(k, s):6.2f}" for k in keys))

P("\nOuter-s drain inside each block: max over the block of (full s - f_s), per flavour (1 = s shell empties completely)")
blocks = {"3d (4s)": ("4s", 10.0, 15.0), "4d (5s)": ("5s", 19.0, 24.0), "4f+5d (6s)": ("6s", 28.0, 40.0), "5f+6d (7s)": ("7s", 44.0, 56.0)}
P("block       " + " ".join(f"{k:>6}" for k in keys))
for b, (s, a, c) in blocks.items():
    row = f"{b:<11} "
    for k in keys:
        z, f, e = D[k]["z"], D[k]["f"][:, lab.index(s)], D[k]["err"]
        w = (z > a + 0.2) & (z < c - 0.2) & (e < 1e-4)
        row += f" {1 - f[w].min() if w.any() else np.nan:5.2f}"
    P(row)

P("\nNon-converged points (err > 1e-4) by z-range:")
for k in keys:
    z, e = D[k]["z"], D[k]["err"]
    bad = z[e > 1e-4]
    P(f"  g={k:>3}: {len(bad)} points" + (f", z in [{bad.min():.1f}, {bad.max():.1f}]" if len(bad) else ""))
open(FN.replace(".npz", "_summary.txt"), "w", encoding="utf-8").write("\n".join(out))
