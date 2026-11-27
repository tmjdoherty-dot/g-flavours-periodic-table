import time
from hf import *
G = Grid(rmin=1e-9, N=1600)
for x, y in ((2, 0), (1, 1), (0, 2)):
    t = time.time()
    sh = core(XE, 2) + [(6, 0, 2), (4, 3, x), (5, 2, y)]
    sh = [s for s in sh if s[2] > 0]
    a = Atom(29, 2, sh, grid=G, qf=hund_qf(sh, 2), lda="fixed")
    o = a.solve(tol=1e-7, maxit=1500)
    i4f = [k for k, s in enumerate(a.shells) if s == (4, 3)]
    print('Ce', x, y, 4 * o['E'], o['it'], o['gmax'], o['conv'], 'eps4f', [4 * o['eps'][k] for k in i4f], f'{time.time()-t:.1f}s', flush=True)
