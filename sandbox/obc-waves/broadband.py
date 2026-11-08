import numpy as np, time, json
from sw1d import *
dy, dt = 10e3, 100.0
yp, ys, sig = 3000e3, 500e3, 100e3
T = 10 * 86400.0
Lext = np.ceil(c * T / dy) * dy
Ly = np.ceil((yp + c * T / 2) / dy) * dy
shape = lambda y: 0.01 * np.exp(-(y - yp)**2 / (2 * sig**2))
gref = Grid(-Lext, Ly, dy); s_ref = southward_packet(gref, shape)
j0 = int(round(Lext / dy)); gm = Grid(0.0, Ly, dy)
s_m = tuple(a[j0:] for a in s_ref)
t, ref, _ = run(gref, s_ref, None, T, dt, ys, out_every=300.0)
alphas = {"flather": c, "tuned_1.6f": admittance(1.6 * f), "half_c": c / 2}
win = np.hanning(len(t))
F = np.fft.rfft(ref * win); w = 2 * np.pi * np.fft.rfftfreq(len(t), t[1] - t[0])
out = {"w_over_f": (w / f).tolist(), "incident": np.abs(F).tolist()}
for name, a in alphas.items():
    _, m, _ = run(gm, s_m, lambda eb, tt, a=a: -a * eb, T, dt, ys, out_every=300.0)
    out[name] = (np.abs(np.fft.rfft((m - ref) * win)) / np.abs(F)).tolist()
    print(name, "done", flush=True)
json.dump(out, open("broadband.json", "w"))
wf = w / f
for target in [1.2, 1.5, 2, 3, 5, 8]:
    i = np.argmin(abs(wf - target))
    print(f"ω/f={wf[i]:.2f}  " + "  ".join(f"{k}: meas {out[k][i]:.3f} th {abs(theory_r(w[i], a)):.3f}" for k, a in alphas.items()))
