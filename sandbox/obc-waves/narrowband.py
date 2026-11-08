import numpy as np, time, json, sys
from sw1d import *

def measure(w0, alphas, dy=10e3, dt=100.0):
    l0 = np.sqrt(w0**2 - f**2) / c
    lam = 2 * np.pi / l0
    sig = 1.5 * lam
    yp = 4 * sig + 500e3
    ys = 500e3
    cg = c**2 * l0 / w0
    T = 1.15 * (yp + 4 * sig + 2 * ys) / cg
    Lext = np.ceil(c * T / dy) * dy
    Ly = np.ceil((yp + 4 * sig + c * T / 2) / dy) * dy
    shape = lambda y: 0.01 * np.exp(-(y - yp)**2 / (2 * sig**2)) * np.cos(l0 * (y - yp))
    gref = Grid(-Lext, Ly, dy)
    s_ref = southward_packet(gref, shape)
    j0 = int(round(Lext / dy))
    gm = Grid(0.0, Ly, dy)
    s_m = (s_ref[0][j0:], s_ref[1][j0:], s_ref[2][j0:])
    _, ref, _ = run(gref, s_ref, None, T, dt, ys)
    out = {}
    for name, a in alphas.items():
        _, m, _ = run(gm, s_m, lambda eb, t, a=a: -a * eb, T, dt, ys)
        out[name] = float(np.sqrt(np.sum((m - ref)**2) / np.sum(ref**2)))
    return out

alphas_for = lambda w: {"flather": c, "tuned_1.6f": admittance(1.6 * f), "dirichlet": 0.0}
res = {}
for r in [1.3, 1.6, 2.0, 3.0, 5.0]:
    t0 = time.time()
    m = measure(r * f, alphas_for(r * f))
    res[r] = m
    th = {k: abs(theory_r(r * f, a)) for k, a in alphas_for(r * f).items()}
    print(f"ω/f = {r}: measured {m}  theory {th}  ({time.time()-t0:.0f} s)", flush=True)
json.dump(res, open("narrowband.json", "w"))
