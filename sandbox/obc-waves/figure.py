# The validation figure from the Python runs (narrowband.py, broadband.py, steady.py, oblique.py,
# rotating_flather.py, oblique_rot.py):  python3 figure.py  ->  obc_1d_validation.png
import numpy as np, json, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sw1d import *
from oblique import theory as theory_oblique
from rotating_flather import alpha_rot, run_rot
nb = json.load(open("narrowband.json")); bb = json.load(open("broadband.json")); st = json.load(open("steady.json"))
rot = json.load(open("rotating_flather.json")); ob = json.load(open("oblique.json")); obr = json.load(open("oblique_rot.json"))

fig, ax = plt.subplots(1, 3, figsize=(19, 5.2), constrained_layout=True)
wf = np.linspace(1.02, 8, 400); w = wf * f
for key, lab, a, col in [("flather", "Flather, α = c", c, "C0"), ("tuned_1.6f", "α = phase speed at 1.6 f", admittance(1.6 * f), "C2"),
                         ("half_c", "α = c/2", c / 2, "C1")]:
    ax[0].plot(wf, np.abs(theory_r(w, a)), color=col, lw=2, label=lab)
    bw = np.array(bb["w_over_f"]); m = (bw >= 2.0) & (bw <= 8)
    ax[0].plot(bw[m], np.convolve(np.array(bb[key]), np.ones(25) / 25, mode="same")[m], color=col, ls="--", lw=1.2)
    if key in nb["1.3"]:
        rs = sorted(float(k) for k in nb); ax[0].plot(rs, [nb[str(r)][key] for r in rs], "o", color=col, ms=7)
ax[0].plot(wf, np.abs((admittance(w) - alpha_rot(w)) / (admittance(w) + alpha_rot(w))), color="k", lw=2,
           label="rotating Flather (new)")
rs = sorted(float(k) for k in rot); ax[0].plot(rs, [rot[str(r)][0] for r in rs], "o", color="k", ms=7)
ax[0].axhline(1, color="C3", lw=2, label="Dirichlet, α = 0")
ax[0].plot([], [], "k--", lw=1, label="runs: broadband pulse"); ax[0].plot([], [], "ko", label="runs: one wave each")
ax[0].axvspan(1.0, 1.27, color="0.9"); ax[0].text(1.3, 0.85, "← longer than a\n   2000 km basin", fontsize=8.5, color="0.4")
ax[0].set(xlabel="wave frequency ω / f", ylabel="reflection coefficient |r|", ylim=(0, 1.05), xlim=(1, 8),
          title="(a) waves leaving at normal incidence")
ax[0].legend(fontsize=8.5, loc="upper right")

th = np.linspace(0, 80, 200)
for r, col in [(1.3, "C0"), (2.0, "C4"), (5.0, "C5")]:
    W = r * f; kap = np.sqrt(W**2 - f**2) / c
    ax[1].plot(th, np.abs(theory_oblique(W, kap * np.sin(np.radians(th)), c)), color=col, lw=2, label=f"Flather, ω = {r} f")
    angs = [0, 30, 60]; ax[1].plot(angs, [ob[f"{r}_{a}"]["flather"][0] for a in angs], "o", color=col, ms=7)
    angs_r = [a for a in (0, 30) if f"{r}_{a}" in obr]
    if angs_r:
        ax[1].plot(angs_r, [obr[f"{r}_{a}"] for a in angs_r], "s", color=col, mfc="none", ms=9, mew=2)
ax[1].plot([], [], "ko", label="runs: Flather"); ax[1].plot([], [], "ks", mfc="none", mew=2, label="runs: rotating Flather")
ax[1].set(xlabel="angle of incidence (°, 0 = straight at the edge)", ylabel="reflection coefficient |r|", ylim=(0, 1),
          xlim=(0, 80), title="(b) waves leaving at an angle")
ax[1].legend(fontsize=8.5, loc="upper left")

y = np.array(st["y"]); gm = Grid(0.0, 3000e3, 10e3)
etaB = np.array(st["truth"]); uB = -(g / f) * np.gradient(etaB, 10e3)
_, (eta_rot, *_ ) = run_rot(gm, (etaB, uB, np.zeros(gm.N + 1)), 5 * 86400.0, 100.0, 500e3,
                            eta_ext=0.1 * np.exp(-(600e3)**2 / (2 * (400e3)**2)) + 0.02, u_ext=uB[0])
ax[2].plot(y / 1e3, 100 * etaB, "k", lw=4, label="truth: steady balanced jet")
ax[2].plot(y / 1e3, 100 * np.array(st["flather_2cm"]), color="C0", lw=2, label="Flather")
ax[2].plot(y / 1e3, 100 * np.array(st["dirichlet_2cm"]), color="C3", lw=2, ls="--", label="Dirichlet (= M5b here)")
ax[2].plot(y / 1e3, 100 * eta_rot, color="C2", lw=2.5, ls=":", label="rotating Flather (new)")
ax[2].set(xlabel="distance from the open edge (km)", ylabel="η (cm)", xlim=(0, 2500),
          title="(c) balanced state, external η 2 cm too high, day 5")
ax[2].legend(fontsize=9)
fig.savefig("obc_1d_validation.png", dpi=150)
print("saved obc_1d_validation.png")
