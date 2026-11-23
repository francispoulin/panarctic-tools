# Scorecard figure for the GitHub issue (claude/issue-draft-mixing-schemes.md): depth error of CATKE, k-ε and NEMO TKE
# against the LES, 4 m columns. Numbers typed in from the comparison tables (22 Nov). python3 scorecard.py -> scorecard.png
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np

# Depth error of each scheme against the LES, % of the LES depth, on the 4 m column grid
groups = [
 ("Convection only (rung 1)", [
   ("Cooling 250 W/m² · N² 1e-5",           [0, 8, -8]),
   ("Cooling 500 W/m² · N² 1e-5",         [6, 6, -6]),
   ("Cooling 1000 W/m² · N² 1e-5",         [4, 0, -12]),
   ("Cooling 1000 W/m² · N² 4e-5",         [0, 8, -8]),
   ("Cooling 250 W/m² · N² 2.5e-6",         [4, 4, -8]),
   ("Arctic: brine of 5 cm/day ice growth",     [9, 0, 0]),
 ]),
 ("Ice–ocean stress, no waves (rungs 2–4)", [
   ("Stress 0.41 N/m², no rotation",     [4, -3, 6]),
   ("Stress 0.41 N/m², rotating (75° N)", [27, 12, -3]),
   ("Arctic: ice drifting at 0.2 m/s",   [40, 27, 0]),
   ("Arctic: drifting ice + brine",      [9, 9, 0]),
   ("Arctic: drifting ice + melt",       [114, 29, 29]),
 ]),
 ("Wind + cooling, no waves (like under ice)", [
   ("Wind 0.31 N/m² + cooling 420 W/m²", [12, 8, -16]),
   ("Wind 0.46 N/m² + cooling 330 W/m²", [21, 9, -13]),
   ("Wind 0.61 N/m² + cooling 210 W/m²", [28, 14, 0]),
   ("Wind 0.70 N/m², no cooling",        [40, 0, 17]),
 ]),
 ("Wind + cooling with waves (open water)", [
   ("Wind 0.31 N/m² + cooling 420 W/m²", [8, 4, 16]),
   ("Wind 0.46 N/m² + cooling 330 W/m²", [16, 4, 20]),
   ("Wind 0.61 N/m² + cooling 210 W/m²", [12, 0, 4]),
   ("Wind 0.70 N/m², no cooling",        [14, -19, -23]),
 ]),
]
schemes = ["CATKE", "k-ε", "NEMO TKE"]

rows, labels, group_rows = [], [], []
for g, items in groups:
    group_rows.append((g, len(rows)))
    for lab, vals in items:
        rows.append(vals); labels.append(lab)
data = np.array(rows, float)

blue, red, mid = "#2a78d6", "#e34948", "#f0efec"
cmap = LinearSegmentedColormap.from_list("div", [blue, mid, red])
norm = TwoSlopeNorm(vmin=-35, vcenter=0, vmax=35)

n = len(rows)
gap = 0.6
ypos, y = [], 0.0
gi = 0
group_y = []
for i in range(n):
    if gi < len(group_rows) and group_rows[gi][1] == i:
        y += gap; group_y.append((group_rows[gi][0], y)); y += 0.55; gi += 1
    ypos.append(y); y += 1.0
H = y

fig, ax = plt.subplots(figsize=(9.2, 0.42 * H + 1.6), dpi=200)
fig.patch.set_facecolor("#fcfcfb"); ax.set_facecolor("#fcfcfb")
ink, muted = "#1f1f1d", "#6b6a66"
for i in range(n):
    for j in range(3):
        v = data[i, j]
        ax.add_patch(plt.Rectangle((j + 0.04, ypos[i] + 0.04), 0.92, 0.92, color=cmap(norm(np.clip(v, -35, 35))),
                                   lw=0))
        mark = "✓ " if abs(v) <= 10 else ""
        ax.text(j + 0.5, ypos[i] + 0.5, (f"{mark}{v:+.0f} %" if v != 0 else "✓ 0 %").replace("-", "−"), ha="center", va="center",
                fontsize=9.5, color=ink, fontweight="bold" if abs(v) > 20 else "normal")
    ax.text(-0.08, ypos[i] + 0.5, labels[i], ha="right", va="center", fontsize=9.5, color=ink)
for g, yy in group_y:
    ax.text(-3.55, yy + 0.25, g, ha="left", va="center", fontsize=10.5, color=ink, fontweight="bold")
for j, s in enumerate(schemes):
    ax.text(j + 0.5, -0.35, s, ha="center", va="center", fontsize=11, color=ink, fontweight="bold")
ax.set_xlim(-3.6, 3.05); ax.set_ylim(H + 0.2, -0.9); ax.axis("off")

# legend bar
cax = fig.add_axes([0.38, 0.045, 0.5, 0.018])
grad = np.linspace(-35, 35, 256)[None, :]
cax.imshow(grad, aspect="auto", cmap=cmap, norm=norm, extent=[-35, 35, 0, 1])
cax.set_yticks([]); cax.set_xticks([-30, -10, 0, 10, 30])
cax.set_xticklabels(["−30 %", "−10 %", "0", "+10 %", "+30 %"], fontsize=8.5, color=muted)
for s in cax.spines.values(): s.set_visible(False)
cax.tick_params(length=0)
fig.text(0.38, 0.075, "Mixed layer too shallow ← → too deep   (✓ = within 10 % of the LES)", fontsize=9, color=muted)
fig.text(0.02, 0.985, "How deep does each scheme mix, compared with the turbulence-resolving simulation (LES)?",
         fontsize=12, color=ink, fontweight="bold", va="top")
fig.text(0.02, 0.957, "Error in mixed-layer depth after 24 h (Arctic cases: 4 days), single column at 4 m. Cooling: heat loss for α = 2×10⁻⁴ K⁻¹;\nwind: stress in N/m² (0.31–0.70 ≈ 13–20 m/s wind); N²: stratification in s⁻² (1e-5 ≈ 0.1 kg/m³ denser per 100 m)",
         fontsize=9.5, color=muted, va="top")
fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.11)
fig.savefig("scorecard.png", facecolor=fig.get_facecolor())
