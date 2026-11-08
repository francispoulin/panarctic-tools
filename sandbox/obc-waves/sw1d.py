"""Linear rotating shallow water, independent of x (along-edge wavenumber k = 0), on a C-grid in y.
South side open (face 0), north side a wall. This is exactly the edge-mean mode of the 2D problem."""
import numpy as np

g, H, f = 9.81, 500.0, 1.4e-4
c = np.sqrt(g * H)

def admittance(w):            # V/η of an outgoing Poincaré wave of frequency w (k = 0): the wave's phase speed
    return c * w / np.sqrt(w**2 - f**2)

def theory_r(w, alpha):       # reflection coefficient (η) of V − Vext = −α (η − ηext)
    Y = admittance(w)
    return (Y - alpha) / (Y + alpha)

class Grid:
    def __init__(self, y0, y1, dy):
        self.N = int(round((y1 - y0) / dy)); self.dy = dy
        self.yc = y0 + (np.arange(self.N) + 0.5) * dy
        self.yf = y0 + np.arange(self.N + 1) * dy

def southward_packet(grid, shape):
    """η, u, v of a purely southward-travelling set of Poincaré waves whose η at t = 0 is shape(yc)."""
    eta0 = shape(grid.yc)
    N, dy = grid.N, grid.dy
    l = 2 * np.pi * np.fft.fftfreq(N, dy)
    eh = np.fft.fft(eta0); eh[0] = 0
    Om = np.sqrt(f**2 + c**2 * l**2)
    with np.errstate(divide="ignore", invalid="ignore"):
        vh = np.where(l != 0, -Om / (H * np.abs(l)) * eh, 0)       # southward: V = −Y η
        om = -np.sign(l) * Om
        uh = np.where(l != 0, 1j * f * vh / om, 0)
    eta = np.real(np.fft.ifft(eh))
    u = np.real(np.fft.ifft(uh))
    vc = np.real(np.fft.ifft(vh * np.exp(-1j * l * dy / 2)))   # values at the south face of each cell
    v = np.append(vc, 0.0)
    return eta, u, v

def run(grid, state0, bc, T, dt, station, ext=lambda t: (0.0, 0.0), out_every=600.0):
    """bc(eta_b, t, extra) -> V at the open face. Returns times and η at the station."""
    eta, u, v = [a.copy() for a in state0]
    N, dy = grid.N, grid.dy
    js = np.argmin(abs(grid.yc - station))
    v[-1] = 0.0

    def rhs(eta, u, v, t):
        v = v.copy()
        if bc is not None:
            v[0] = bc(eta[0], t) / H
        vbar = 0.5 * (v[:-1] + v[1:])
        ubar = np.zeros(N + 1); ubar[1:-1] = 0.5 * (u[:-1] + u[1:])
        deta = -H * (v[1:] - v[:-1]) / dy
        du = f * vbar
        dv = np.zeros(N + 1)
        dv[1:-1] = -f * ubar[1:-1] - g * (eta[1:] - eta[:-1]) / dy
        return deta, du, dv, v

    nsteps = int(round(T / dt)); every = int(round(out_every / dt))
    ts, rec = [0.0], [eta[js]]
    t = 0.0
    for n in range(nsteps):
        k1 = rhs(eta, u, v, t)
        k2 = rhs(eta + dt/2*k1[0], u + dt/2*k1[1], v + dt/2*k1[2], t + dt/2)
        k3 = rhs(eta + dt/2*k2[0], u + dt/2*k2[1], v + dt/2*k2[2], t + dt/2)
        k4 = rhs(eta + dt*k3[0], u + dt*k3[1], v + dt*k3[2], t + dt)
        eta = eta + dt/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])
        u   = u   + dt/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
        v   = v   + dt/6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])
        t += dt
        if (n + 1) % every == 0:
            ts.append(t); rec.append(eta[js])
    if bc is not None:
        v[0] = bc(eta[0], t) / H
    return np.array(ts), np.array(rec), (eta, u, v)
