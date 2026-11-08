import numpy as np, json
from sw1d import *
dy, dt, Ly, T = 10e3, 100.0, 3000e3, 5 * 86400.0
gm = Grid(0.0, Ly, dy); y = gm.yc
yj, wj, A = 600e3, 400e3, 0.1
etaB = A * np.exp(-(y - yj)**2 / (2 * wj**2))
uB = -(g / f) * np.gradient(etaB, dy)           # geostrophic along-edge jet; v = 0
eB0 = A * np.exp(-(0 - yj)**2 / (2 * wj**2))    # balanced η at the open face
s0 = (etaB, uB, np.zeros(gm.N + 1))
out = {"y": y.tolist(), "truth": etaB.tolist()}
for delta in [0.0, 0.02]:
    for name, a in {"flather": c, "dirichlet": 0.0}.items():
        _, _, (eta, u, v) = run(gm, s0, lambda eb, t, a=a: -a * (eb - (eB0 + delta)), T, dt, 500e3)
        err = np.sqrt(np.mean((eta - etaB)**2))
        jet = H * np.sum(u - uB) * dy / 1e6
        print(f"datum {100*delta:.0f} cm, {name:9s}: rms η error {100*err:.3f} cm, spurious along-edge transport {jet:+.3f} Sv")
        out[f"{name}_{int(100*delta)}cm"] = eta.tolist()
print("gHδ/f for δ = 2 cm:", g * H * 0.02 / f / 1e6, "Sv")
json.dump(out, open("steady.json", "w"))
