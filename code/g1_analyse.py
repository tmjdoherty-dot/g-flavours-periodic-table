"""Analyse results/g1_scan.npz: filling order vs Madelung for each g (lam = 1/g)."""
import numpy as np

import sys
d = np.load(sys.argv[1] if len(sys.argv) > 1 else "results/g1_scan_run2.npz", allow_pickle=True)
LN = "spdfg"
labels = [tuple(x) for x in d["labels"]]
names = [f"{n}{LN[l]}" for n, l in labels]
cap = np.array([2 * l + 1 for n, l in labels], float)
ORDER = sorted(range(len(labels)), key=lambda i: (labels[i][0] + labels[i][1], labels[i][0]))


def madelung(z):
    f = np.zeros(len(labels))
    left = z
    for i in ORDER:
        take = min(cap[i], left)
        f[i] = take
        left -= take
        if left <= 0:
            break
    return f


# Madelung onset (per flavour) = electrons before the shell starts
cum, mad_on = 0.0, {}
for i in ORDER:
    mad_on[names[i]] = cum
    cum += cap[i]

SHELLS = ["2p", "3s", "3p", "4s", "3d", "4p", "5s", "4d", "5p", "6s", "4f", "5d", "6p", "7s", "5f", "6d", "7p", "8s", "5g"]
EXC = {24, 29, 41, 42, 44, 45, 46, 47, 57, 58, 64, 78, 79, 89, 90, 91, 92, 93, 96, 103}

out = []
P = lambda s="": (print(s), out.append(s))
P("Onsets (first z with f > 0.02 per flavour) vs Madelung; z = Z/g. Real elements: Z = g z.")
keys = ["inf", "4", "3", "2", "2c", "1"]
P("shell  Mad   " + "  ".join(f"{k:>6}" for k in keys))
data = {}
for k in keys:
    m = d["key"] == k
    o = np.argsort(d["z"][m])
    data[k] = dict(z=d["z"][m][o], f=d["f"][m][o], err=d["err"][m][o], mu=d["mu"][m][o], eps=d["eps"][m][o])
onsets = {k: {} for k in keys}
for s in SHELLS:
    i = names.index(s)
    row = f"{s:>5} {mad_on[s]:5.2f} "
    for k in keys:
        z, f = data[k]["z"], data[k]["f"]
        hit = np.nonzero(f[:, i] > 0.02)[0]
        on = z[hit[0]] if len(hit) else np.nan
        onsets[k][s] = on
        row += f"  {on:6.2f}"
    P(row)

P("\nFilling order by onset (ties broken by Madelung order):")
for k in keys:
    seq = sorted([s for s in SHELLS if np.isfinite(onsets[k][s])], key=lambda s: (onsets[k][s], mad_on[s]))
    P(f"  g={k:>3}: " + " ".join(seq))
P("  Madelung: " + " ".join(sorted(SHELLS, key=lambda s: mad_on[s])))

P("\nDeviation from Madelung: share of the z-range [0.25, 60] where |f - f_Madelung|_1 > 0.1 per flavour")
for k in keys:
    z, f = data[k]["z"], data[k]["f"]
    dev = np.array([np.abs(fi - madelung(zi)).sum() for zi, fi in zip(z, f)])
    bad = dev > 0.1
    P(f"  g={k:>3}: {bad.mean():.3f}  (unconverged points err>1e-5: {(data[k]['err'] > 1e-5).sum()})")
    data[k]["dev"] = dev

P("\nReal chemistry check (g = 2 rows at integer Z): non-Madelung if total |f - f_Mad|_1 > 1 electron (Cr-type 4s->3d shift = 2)")
for k in ("2", "2c"):
    z, f = data[k]["z"], data[k]["f"]
    Zs = np.round(2 * z).astype(int)
    sel = np.isclose(2 * z, Zs)
    nm = {Z for Z, zi, fi in zip(Zs[sel], z[sel], f[sel]) if np.abs(fi - madelung(zi)).sum() * 2 > 1.0}
    P(f"  {k}: LDA non-Madelung Z = {sorted(nm)}")
    P(f"      experiment exceptions caught {len(nm & EXC)}/{len(EXC)}; extra (LDA only) {len(nm - EXC)}")

open("results/g1_summary.txt", "w", encoding="utf-8").write("\n".join(out))
np.savez("results/g1_onsets.npz", **{k: np.array([onsets[k][s] for s in SHELLS]) for k in keys}, shells=np.array(SHELLS))
