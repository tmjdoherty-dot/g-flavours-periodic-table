"""G2 part C: lanthanide row, Hund HF: [Xe] 6s^g 4f^x 5d^y, y in {0,1,2}. Ground y per element.
La/Ac step = the first element (delta = 1) takes 5d (y = 1) before 4f.
Real g = 2 (experiment): La 5d1, Ce 4f1 5d1, Pr..Eu 4f^n, Gd 4f7 5d1, Tb..Yb 4f^n, Lu 4f14 5d1."""
import json, sys
from collections import defaultdict
R = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "results/g2/C.json"))
E = defaultdict(dict)
for r in R:
    if "err" in r or not r["conv"]:
        continue
    E[(r["g"], r["delta"], r["x"], r["y"])][r["N"]] = r["E"]
ext = defaultdict(dict)
for (g, d, x, y), v in E.items():
    if 1600 in v and 2400 in v:
        ext[(g, d)][y] = v[2400] + (v[2400] - v[1600]) / 1.25
out = []
P = lambda s="": (print(s), out.append(s))
P(f"rows {len(R)}, failed/unconverged {sum(1 for r in R if 'err' in r or not r['conv'])}")
EXP = [1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]  # La..Yb ground 5d count (delta 1..14)
for g in sorted({k[0] for k in ext}):
    ys = []
    for d in range(1, 7 * g + 1):
        e = ext.get((g, d), {})
        ys.append(str(min(e, key=e.get)) if e else "?")
    first = ext.get((g, 1), {})
    gap = (first.get(0, float("nan")) - first.get(1, float("nan"))) * g * g * 27.2114
    nd = sum(1 for c in ys if c not in ("0", "?"))
    P(f"g={g}: 5d count of ground per element: {' '.join(ys)}")
    P(f"      La-type step (delta=1 takes 5d): {'yes' if ys[0] == '1' else 'no'}  "
      f"E(4f1)-E(5d1) = {gap:+.2f} eV;  elements with 5d occupied: {nd}/{7*g}")
    if g == 2:
        P(f"      experiment: {' '.join(map(str, EXP))};  agree {sum(1 for a, b in zip(ys, EXP) if a == str(b))}/14")
open("results/g2/C_summary.txt", "w", encoding="utf-8").write("\n".join(out))
