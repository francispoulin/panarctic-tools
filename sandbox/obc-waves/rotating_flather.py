"""A rotation-aware Flather condition for the edge-mean mode (k = 0), tested like Flather.

An outgoing Poincaré wave has V/η = cₚ = c/√(1 − f²/ω²); Flather uses c. To second order in f²/ω²,
cₚ ≈ c/(1 − f²/2ω²), i.e.  V (1 − f²/2ω²) = −c η.  For k = 0 waves ∂u/∂t = f v, so (f²/ω²) V = −f H ∫ u dt, and the
condition becomes local in time:
    V = Vᵉˣᵗ − c (ηᵇ − ηᵉˣᵗ) + S/2,      dS/dt = −f H (uᵇ − uᵉˣᵗ),   S(0) = 0.
S only sees the along-edge velocity, so it is blind to a datum in ηᵉˣᵗ (the datum still enters through −c(ηᵇ − ηᵉˣᵗ)).
Effective weight α = c/(1 − f²/2ω²), against the exact cₚ."""
import numpy as np, json
from sw1d import g, H, f, c, Grid, southward_packet, admittance

def alpha_rot(w, a=0.5):
    return c / (1 - a * f**2 / w**2)

def run_rot(grid, s0, T, dt, ys, a=0.5, eta_ext=0.0, u_ext=0.0, out_every=600.0):
    eta, u, v = [x.copy() for x in s0]
    S = 0.0
    N, dy = grid.N, grid.dy
    js = np.argmin(abs(grid.yc - ys))
    def rhs(eta, u, v, S):
        v = v.copy()
        v[0] = (-c * (eta[0] - eta_ext) + a * S) / H
        vbar = 0.5 * (v[:-1] + v[1:])
        ubar = np.zeros(N + 1); ubar[1:-1] = 0.5 * (u[:-1] + u[1:])
        dv = np.zeros(N + 1)
        dv[1:-1] = -f * ubar[1:-1] - g * (eta[1:] - eta[:-1]) / dy
        return -H * (v[1:] - v[:-1]) / dy, f * vbar, dv, -f * H * (u[0] - u_ext)
    every = int(round(out_every / dt)); rec = [eta[js]]
    for n in range(int(round(T / dt))):
        st = (eta, u, v, S)
        k1 = rhs(*st)
        k2 = rhs(*(x + dt / 2 * y for x, y in zip(st, k1)))
        k3 = rhs(*(x + dt / 2 * y for x, y in zip(st, k2)))
        k4 = rhs(*(x + dt * y for x, y in zip(st, k3)))
        eta, u, v, S = (x + dt / 6 * (y1 + 2 * y2 + 2 * y3 + y4) for x, y1, y2, y3, y4 in zip(st, k1, k2, k3, k4))
        v[-1] = 0
        if (n + 1) % every == 0:
            rec.append(eta[js])
    return np.array(rec), (eta, u, v, S)

if __name__ == "__main__":
    from sw1d import run
    dy, dt, ys = 10e3, 100.0, 500e3
    res = {}
    for r in (1.3, 1.6, 2.0, 3.0, 5.0):
        w0 = r * f; l0 = np.sqrt(w0**2 - f**2) / c; sig = 1.5 * 2 * np.pi / l0; yp = 4 * sig + 500e3
        cg = c**2 * l0 / w0; T = 1.15 * (yp + 4 * sig + 2 * ys) / cg
        Lext = np.ceil(c * T / dy) * dy; Ly = np.ceil((yp + 4 * sig + c * T / 2) / dy) * dy
        shape = lambda y: 0.01 * np.exp(-(y - yp)**2 / (2 * sig**2)) * np.cos(l0 * (y - yp))
        gref = Grid(-Lext, Ly, dy); s_ref = southward_packet(gref, shape)
        j0 = int(round(Lext / dy)); gm = Grid(0.0, Ly, dy); s_m = tuple(x[j0:] for x in s_ref)
        _, ref, _ = run(gref, s_ref, None, T, dt, ys)
        m, _ = run_rot(gm, s_m, T, dt, ys)
        meas = float(np.sqrt(np.sum((m - ref)**2) / np.sum(ref**2)))
        Y = admittance(w0); th = abs((Y - alpha_rot(w0)) / (Y + alpha_rot(w0)))
        res[r] = (meas, float(th))
        print(f"ω/f = {r}: rotating Flather measured {meas:.3f}, theory {th:.3f}   (Flather theory {abs((Y-c)/(Y+c)):.3f})",
              flush=True)
    # rung 0: steady jet with a 2 cm datum
    gm = Grid(0.0, 3000e3, dy); y = gm.yc
    etaB = 0.1 * np.exp(-(y - 600e3)**2 / (2 * (400e3)**2)); uB = -(g / f) * np.gradient(etaB, dy)
    for delta in (0.0, 0.02):
        _, (eta, u, v, S) = run_rot(gm, (etaB, uB, np.zeros(gm.N + 1)), 5 * 86400.0, dt, ys,
                                    eta_ext=0.1 * np.exp(-(600e3)**2 / (2 * (400e3)**2)) + delta, u_ext=uB[0])
        print(f"steady jet, datum {100*delta:.0f} cm: rms η error {100*np.sqrt(np.mean((eta-etaB)**2)):.3f} cm, "
              f"spurious transport {H*np.sum(u-uB)*dy/1e6:+.3f} Sv")
    json.dump(res, open("rotating_flather.json", "w"))
