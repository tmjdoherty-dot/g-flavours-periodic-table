"""Validate atom.py against NIST LDA reference (VWN; PZ81 used here, differences ~1e-4 Ha)."""
import time, numpy as np
from atom import Grid, scf, config_str
grid = Grid()
ref = {10: {(1,0): -30.305855, (2,0): -1.322610, (2,1): -0.498034},
       18: {(1,0): -113.800134, (2,0): -10.794172, (2,1): -8.443439, (3,0): -0.883384, (3,1): -0.382330}}
for Z in (10, 18):
    t = time.time()
    out = scf(Z/2, 0.5, grid, g_corr=2, kT=1e-4)
    print(f"Z={Z} it={out['it']} err={out['err']:.1e} t={time.time()-t:.2f}s")
    for k, e in ref[Z].items():
        print("  ", k, f"{4*out['occ'][k][1]:.5f}  NIST {e:.5f}")
