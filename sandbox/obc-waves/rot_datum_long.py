import numpy as np
from sw1d import g, H, f, c, Grid, run
from rotating_flather import run_rot
dy, dt = 10e3, 100.0
gm = Grid(0.0, 3000e3, dy); y = gm.yc
etaB = 0.1 * np.exp(-(y - 600e3)**2 / (2 * (400e3)**2)); uB = -(g / f) * np.gradient(etaB, dy)
eext = 0.1 * np.exp(-(600e3)**2 / (2 * (400e3)**2)) + 0.02
for days in (0.25, 1, 5, 30):
    _, (eta, u, v, S) = run_rot(gm, (etaB, uB, np.zeros(gm.N + 1)), days * 86400.0, dt, 500e3, eta_ext=eext, u_ext=uB[0])
    _, _, (e2, u2, v2) = run(gm, (etaB, uB, np.zeros(gm.N + 1)), lambda eb, t: -c * (eb - eext), days * 86400.0, dt, 500e3)
    print(f"day {days:5}: rotating Flather rms err {100*np.sqrt(np.mean((eta-etaB)**2)):.3f} cm, S/2c = {100*S/2/c:+.3f} cm | "
          f"Flather rms err {100*np.sqrt(np.mean((e2-etaB)**2)):.3f} cm", flush=True)
