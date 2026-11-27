"""Validate hf.py at g = 2 against numerical HF reference energies (Clementi-Roetti / Bunge tables)."""
import time
from hf import Atom, core, AR
REF = {"Be": (4, [(1,0,2),(2,0,2)], -14.573023), "Ne": (10, core([(1,0),(2,0),(2,1)],2), -128.547098),
       "Na": (11, core([(1,0),(2,0),(2,1)],2)+[(3,0,1)], -161.858912),
       "Ar": (18, core(AR,2), -526.817513), "K": (19, core(AR,2)+[(4,0,1)], -599.164787),
       "Cu": (29, core(AR,2)+[(3,2,10),(4,0,1)], -1638.963742), "Zn": (30, core(AR,2)+[(3,2,10),(4,0,2)], -1777.848116)}
for name, (Z, sh, ref) in REF.items():
    t = time.time()
    a = Atom(Z/2, 2, sh)
    o = a.solve(tol=1e-6)
    print(f"{name:3s} E={4*o['E']:.6f} ref={ref:.6f} diff={4*o['E']-ref:+.5f} it={o['it']} conv={o['conv']} {time.time()-t:.1f}s", flush=True)
