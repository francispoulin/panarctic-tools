"""Rotating Flather (V = −cη + S/2, dS/dt = −f H uᵇ) for waves arriving at an angle (k ≠ 0)."""
import numpy as np, json
from sw1d import g, H, f, c, Grid
from oblique import packet, theory

def run_rot(grid, k, s0, T, dt, ys, a=0.5):
    eta, u, v = [x.astype(complex).copy() for x in s0]; S = 0j
    N, dy = grid.N, grid.dy; js = np.argmin(abs(grid.yc - ys))
    def rhs(eta, u, v, S):
        v = v.copy(); v[0] = (-c * eta[0] + a * S) / H
        vbar = 0.5 * (v[:-1] + v[1:]); ubar = np.zeros(N + 1, complex); ubar[1:-1] = 0.5 * (u[:-1] + u[1:])
        dv = np.zeros(N + 1, complex); dv[1:-1] = -f * ubar[1:-1] - g * (eta[1:] - eta[:-1]) / dy
        return -H * (1j * k * u + (v[1:] - v[:-1]) / dy), f * vbar - 1j * k * g * eta, dv, -f * H * u[0]
    rec = [eta[js]]
    for n in range(int(round(T / dt))):
        st = (eta, u, v, S)
        k1 = rhs(*st); k2 = rhs(*(x + dt/2*y for x, y in zip(st, k1)))
        k3 = rhs(*(x + dt/2*y for x, y in zip(st, k2))); k4 = rhs(*(x + dt*y for x, y in zip(st, k3)))
        eta, u, v, S = (x + dt/6*(a1 + 2*a2 + 2*a3 + a4) for x, a1, a2, a3, a4 in zip(st, k1, k2, k3, k4))
        v[-1] = 0
        if (n + 1) % 6 == 0: rec.append(eta[js])
    return np.array(rec)

from oblique import run as run_ref
res = {}
dy, dt, ys = 10e3, 100.0, 500e3
for r in (1.3, 2.0):
    for th in (0, 30):
        w = r * f; kappa = np.sqrt(w**2 - f**2) / c
        k = kappa * np.sin(np.radians(th)); l = kappa * np.cos(np.radians(th))
        sig = 1.0 * 2 * np.pi / l; yp = 4 * sig + 500e3; vg = c**2 * l / w
        T = np.ceil(1.15 * (yp + 4 * sig + 2 * ys) / vg / (6 * dt)) * 6 * dt
        Lext = np.ceil(c * T / dy) * dy; Ly = np.ceil((yp + 4 * sig + c * T / 2) / dy) * dy
        gref = Grid(-Lext, Ly, dy); s = packet(gref, w, k, yp, sig)
        j0 = int(round(Lext / dy)); gm = Grid(0.0, Ly, dy)
        ref = run_ref(gref, k, s, 0.0, T, dt, ys, open_south=False)
        m = run_rot(gm, k, tuple(x[j0:] for x in s), T, dt, ys)
        meas = float(np.sqrt(np.sum(abs(m - ref)**2) / np.sum(abs(ref)**2)))
        res[f"{r}_{th}"] = meas
        print(f"ω/f = {r}, angle {th}°: rotating Flather {meas:.3f}  (Flather theory {abs(theory(w, k, c)):.3f})", flush=True)
json.dump(res, open("oblique_rot.json", "w"))
