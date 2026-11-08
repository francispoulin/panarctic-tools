"""Rung 3: waves that reach the open edge at an angle. Fields vary as exp(ikx) along the edge (k ≠ 0), so the
problem stays one-dimensional in y but complex. Measured reflection against theory for V − Vext = −α(η − ηext):
    outgoing (southward)  V = −Yₛ η,  Yₛ = gH (ω|l| + i f k) / (ω² − f²)
    reflected (northward) V = +Yₙ η,  Yₙ = gH (ω|l| − i f k) / (ω² − f²)
    r = (Yₛ − α) / (Yₙ + α),  with  ω² = f² + c²(k² + l²)."""
import numpy as np, json
from sw1d import g, H, f, c, Grid

def theory(w, k, alpha):
    l = np.sqrt((w**2 - f**2) / c**2 - k**2)
    Ys = g * H * (w * l + 1j * f * k) / (w**2 - f**2)
    Yn = g * H * (w * l - 1j * f * k) / (w**2 - f**2)
    return (Ys - alpha) / (Yn + alpha)

def run(grid, k, s0, alpha, T, dt, ys, open_south=True):
    eta, u, v = [a.astype(complex).copy() for a in s0]
    N, dy = grid.N, grid.dy
    js = np.argmin(abs(grid.yc - ys))
    def rhs(eta, u, v):
        v = v.copy()
        if open_south:
            v[0] = -alpha * eta[0] / H
        vbar = 0.5 * (v[:-1] + v[1:])
        ubar = np.zeros(N + 1, complex); ubar[1:-1] = 0.5 * (u[:-1] + u[1:])
        deta = -H * (1j * k * u + (v[1:] - v[:-1]) / dy)
        du = f * vbar - 1j * k * g * eta
        dv = np.zeros(N + 1, complex)
        dv[1:-1] = -f * ubar[1:-1] - g * (eta[1:] - eta[:-1]) / dy
        return deta, du, dv
    rec = [eta[js]]
    for n in range(int(round(T / dt))):
        k1 = rhs(eta, u, v)
        k2 = rhs(*(a + dt / 2 * b for a, b in zip((eta, u, v), k1)))
        k3 = rhs(*(a + dt / 2 * b for a, b in zip((eta, u, v), k2)))
        k4 = rhs(*(a + dt * b for a, b in zip((eta, u, v), k3)))
        eta, u, v = (a + dt / 6 * (b1 + 2 * b2 + 2 * b3 + b4) for a, b1, b2, b3, b4 in zip((eta, u, v), k1, k2, k3, k4))
        v[-1] = 0
        if (n + 1) % 6 == 0:
            rec.append(eta[js])
    return np.array(rec)

def packet(grid, w, k, yp, sig, A=0.01):
    l = -np.sqrt((w**2 - f**2) / c**2 - k**2)                 # southward
    Ys = g * H * (w * -l + 1j * f * k) / (w**2 - f**2)
    env = lambda y: A * np.exp(-(y - yp)**2 / (2 * sig**2)) * np.exp(1j * l * (y - yp))
    eta = env(grid.yc)
    vv = -(Ys / H) * env(grid.yf); vv[-1] = 0
    vc = -(Ys / H) * eta
    u = (1j * f * vc + g * k * eta) / w
    return eta, u, vv

def measure(w, theta_deg, alphas, dy=10e3, dt=100.0, ys=500e3):
    kappa = np.sqrt(w**2 - f**2) / c
    k = kappa * np.sin(np.radians(theta_deg)); l = kappa * np.cos(np.radians(theta_deg))
    lam = 2 * np.pi / l; sig = 1.0 * lam; yp = 4 * sig + 500e3
    vg = c**2 * l / w                                          # group velocity towards the edge
    T = 1.15 * (yp + 4 * sig + 2 * ys) / vg
    T = np.ceil(T / (6 * dt)) * 6 * dt
    Lext = np.ceil(c * T / dy) * dy; Ly = np.ceil((yp + 4 * sig + c * T / 2) / dy) * dy
    gref = Grid(-Lext, Ly, dy); s = packet(gref, w, k, yp, sig)
    j0 = int(round(Lext / dy)); gm = Grid(0.0, Ly, dy)
    ref = run(gref, k, s, 0.0, T, dt, ys, open_south=False)
    out = {}
    for name, a in alphas.items():
        m = run(gm, k, tuple(x[j0:] for x in s), a, T, dt, ys)
        out[name] = (float(np.sqrt(np.sum(abs(m - ref)**2) / np.sum(abs(ref)**2))), float(abs(theory(w, k, a))))
    return out

if __name__ == "__main__":
    res = {}
    for r in (1.3, 2.0, 5.0):
        for th in (0, 30, 60):
            o = measure(r * f, th, {"flather": c, "dirichlet": 0.0})
            res[f"{r}_{th}"] = o
            print(f"ω/f = {r}, angle {th:2d}°:  Flather measured {o['flather'][0]:.3f}, theory {o['flather'][1]:.3f};  "
                  f"Dirichlet {o['dirichlet'][0]:.3f} ({o['dirichlet'][1]:.3f})", flush=True)
    json.dump(res, open("oblique.json", "w"))
