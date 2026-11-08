# Check the analytic (envelope × carrier) southward packet used in the Oceananigans script against the spectral one.
import numpy as np
from sw1d import *
def analytic(grid, w0, yp, sig, A=0.01):
    l0 = np.sqrt(w0**2 - f**2) / c; Y = w0 / l0
    env = lambda y: A * np.exp(-(y - yp)**2 / (2 * sig**2))
    th = lambda y: l0 * (y - yp)
    eta = env(grid.yc) * np.cos(th(grid.yc))
    u = -(f * Y / (w0 * H)) * env(grid.yc) * np.sin(th(grid.yc))
    v = -(Y / H) * env(grid.yf) * np.cos(th(grid.yf)); v[-1] = 0
    return eta, u, v
dy, dt = 10e3, 100.0
for r in [1.3, 2.0, 5.0]:
    w0 = r * f; l0 = np.sqrt(w0**2 - f**2) / c; lam = 2*np.pi/l0; sig = 1.5*lam; yp = 4*sig + 500e3; ys = 500e3
    cg = c**2 * l0 / w0; T = 1.15 * (yp + 4*sig + 2*ys) / cg
    Lext = np.ceil(c*T/dy)*dy; Ly = np.ceil((yp + 4*sig + c*T/2)/dy)*dy
    gref = Grid(-Lext, Ly, dy); s_ref = analytic(gref, w0, yp, sig)
    j0 = int(round(Lext/dy)); gm = Grid(0.0, Ly, dy); s_m = tuple(a[j0:] for a in s_ref)
    _, ref, (e, _, _) = run(gref, s_ref, None, T, dt, ys)
    _, m, _ = run(gm, s_m, lambda eb, t: -c*eb, T, dt, ys)
    north = np.sqrt(np.sum(e[gref.yc > yp + 5*sig]**2) / np.sum(s_ref[0]**2))
    print(f"ω/f={r}: Flather r = {np.sqrt(np.sum((m-ref)**2)/np.sum(ref**2)):.3f} (theory {abs(theory_r(w0, c)):.3f}); "
          f"energy leaking north from the analytic start {north:.3f}")
