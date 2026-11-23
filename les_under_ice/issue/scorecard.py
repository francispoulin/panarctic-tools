# Scorecard figures for the Oceananigans issue (claude/issue-draft-catke-under-ice.md), from scorecard_numbers.csv
# (written by scorecard_numbers.jl on pelee, 23 Nov: every depth with the current, interpolated measure).
#   scorecard.png        CATKE, k-ε and NEMO's TKE against the LES, 4 m columns, Δt = 1 min
#   scorecard_catke.png  default against adapted CATKE (rotation limit C^f = 0.1, Cˢ = 0.991, half the surface TKE flux),
#                        4 m and GLORYS's 50 levels, Δt = 10 min
# python3 scorecard.py
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np

errors = {}
for r in csv.DictReader(open("scorecard_numbers.csv")):
    errors[(r["case"], r["column"])] = 100 * float(r["error"])

conv = lambda N2, Jb: f"convection_N2_{N2}_Jb_{Jb}"
groups = [
 ("Convection only", [
   (conv("1.0e-5", "1.2e-7"), "Cooling 250 W/m² · N² 1e-5"),
   (conv("1.0e-5", "2.4e-7"), "Cooling 500 W/m² · N² 1e-5"),
   (conv("1.0e-5", "4.8e-7"), "Cooling 1000 W/m² · N² 1e-5"),
   (conv("4.0e-5", "4.8e-7"), "Cooling 1000 W/m² · N² 4e-5"),
   (conv("2.5e-6", "1.2e-7"), "Cooling 250 W/m² · N² 2.5e-6"),
   ("arctic_brine_convection", "Arctic: salt from 5 cm/day ice growth"),
 ]),
 ("Stress and drifting ice, no waves", [
   ("kato_phillips",         "Stress 0.41 N/m², no rotation"),
   ("rotating_stress",       "Stress 0.41 N/m², rotating (75° N)"),
   ("arctic_ice_drag",       "Arctic: ice drifting at 0.2 m/s"),
   ("arctic_ice_drag_brine", "Arctic: drifting ice + freezing"),
   ("arctic_ice_drag_melt",  "Arctic: drifting ice + melt"),
 ]),
 ("Wind + cooling, no waves (like under ice)", [
   ("weak_wind_strong_cooling", "Wind 0.31 N/m² + cooling 420 W/m²"),
   ("med_wind_med_cooling",     "Wind 0.46 N/m² + cooling 330 W/m²"),
   ("strong_wind_weak_cooling", "Wind 0.61 N/m² + cooling 210 W/m²"),
   ("strong_wind",              "Wind 0.70 N/m², no cooling"),
 ]),
 ("Wind + cooling with waves (open water)", [
   ("weak_wind_strong_cooling_stokes", "Wind 0.31 N/m² + cooling 420 W/m²"),
   ("med_wind_med_cooling_stokes",     "Wind 0.46 N/m² + cooling 330 W/m²"),
   ("strong_wind_weak_cooling_stokes", "Wind 0.61 N/m² + cooling 210 W/m²"),
   ("strong_wind_stokes",              "Wind 0.70 N/m², no cooling"),
 ]),
]
targets = {"rotating_stress", "arctic_ice_drag", "arctic_ice_drag_brine", "arctic_ice_drag_melt",
           "weak_wind_strong_cooling", "med_wind_med_cooling", "strong_wind_weak_cooling", "strong_wind"}

blue, red, mid = "#2a78d6", "#e34948", "#f0efec"
cmap = LinearSegmentedColormap.from_list("div", [blue, mid, red])
norm = TwoSlopeNorm(vmin=-35, vcenter=0, vmax=35)
ink, muted, paper = "#1f1f1d", "#6b6a66", "#fcfcfb"


def draw(columns, headers, title, subtitle, filename, mark_targets=False, separator=None):
    # keep only rows with every column present
    kept = [(g, [(c, lab) for c, lab in items if all(np.isfinite(errors.get((c, col), np.nan)) for col in columns)])
            for g, items in groups]
    kept = [(g, items) for g, items in kept if items]
    ncol = len(columns)
    rows, labels, cases, group_y, ypos, y = [], [], [], [], [], 0.0
    for g, items in kept:
        y += 0.6; group_y.append((g, y)); y += 0.55
        for c, lab in items:
            rows.append([errors[(c, col)] for col in columns]); labels.append(lab); cases.append(c)
            ypos.append(y); y += 1.0
    H, data = y, np.array(rows)
    label_width = 3.6
    fig, ax = plt.subplots(figsize=(5.6 + 1.2 * ncol, 0.42 * H + 2.0), dpi=200)
    fig.patch.set_facecolor(paper); ax.set_facecolor(paper)
    for i in range(len(rows)):
        for j in range(ncol):
            v = data[i, j]
            ax.add_patch(plt.Rectangle((j + 0.04, ypos[i] + 0.04), 0.92, 0.92,
                                       color=cmap(norm(np.clip(v, -35, 35))), lw=0))
            r = round(v)
            mark = "✓ " if abs(v) <= 10 else ""
            text = (f"{mark}{r:+d} %" if r != 0 else "✓ 0 %").replace("-", "−")
            ax.text(j + 0.5, ypos[i] + 0.5, text, ha="center", va="center", fontsize=9.5, color=ink,
                    fontweight="bold" if abs(v) > 20 else "normal")
        lab = labels[i] + ("   ◆" if mark_targets and cases[i] in targets else "")
        ax.text(-0.08, ypos[i] + 0.5, lab, ha="right", va="center", fontsize=9.5, color=ink)
    for g, yy in group_y:
        ax.text(-label_width + 0.05, yy + 0.25, g, ha="left", va="center", fontsize=10.5, color=ink, fontweight="bold")
    for j, h in enumerate(headers):
        ax.text(j + 0.5, -0.55, h, ha="center", va="center", fontsize=10, color=ink, fontweight="bold",
                linespacing=1.15)
    if separator is not None:
        ax.plot([separator, separator], [-1.1, H], color=muted, lw=1.2)
    ax.set_xlim(-label_width, ncol + 0.05); ax.set_ylim(H + 0.2, -1.3); ax.axis("off")

    cax = fig.add_axes([0.40, 0.04, 0.48, 0.016])
    grad = np.linspace(-35, 35, 256)[None, :]
    cax.imshow(grad, aspect="auto", cmap=cmap, norm=norm, extent=[-35, 35, 0, 1])
    cax.set_yticks([]); cax.set_xticks([-30, -10, 0, 10, 30])
    cax.set_xticklabels(["−30 %", "−10 %", "0", "+10 %", "+30 %"], fontsize=8.5, color=muted)
    for s in cax.spines.values(): s.set_visible(False)
    cax.tick_params(length=0)
    legend = "Mixed layer too shallow ← → too deep   (✓ = within 10 % of the LES)"
    if mark_targets:
        legend += "\n◆ = target (default CATKE wrong); others are guards (must not get worse)"
    fig.text(0.40, 0.068, legend, fontsize=8.5, color=muted, va="bottom")
    fig.text(0.02, 0.988, title, fontsize=12, color=ink, fontweight="bold", va="top")
    fig.text(0.02, 0.962, subtitle, fontsize=9, color=muted, va="top")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.885, bottom=0.12)
    fig.savefig(filename, facecolor=fig.get_facecolor())


units = ("Error in mixed-layer depth after 1 day (Arctic cases: 4 days). Cooling: heat loss for α = 2×10⁻⁴ K⁻¹; wind:\n"
         "stress in N/m² (0.31–0.70 ≈ 13–20 m/s winds); N²: stratification in s⁻² (1e-5 ≈ 0.1 kg/m³ denser per 100 m)")

draw(["catke_dz4_dt1", "k_epsilon_dz4_dt1", "nemo_dz4_dt1"], ["CATKE", "k-ε", "NEMO TKE"],
     "How deep does each scheme mix, compared with the turbulence-resolving simulation (LES)?",
     units + "\nSingle columns at 4 m, Δt = 1 min. NEMO's wave terms off without waves, on with waves.",
     "scorecard.png")

draw(["catke_stock_dz4_dt10", "catke_adapted_dz4_dt10", "catke_stock_glorys50_dt10", "catke_adapted_glorys50_dt10"],
     ["Default\n4 m", "Adapted\n4 m", "Default\nGLORYS levels", "Adapted\nGLORYS levels"],
     "Adapting CATKE to no waves and rotation: default against adapted",
     units + "\nΔt = 10 min. Adapted: rotation limit ℓ ≤ 0.1 √e/|f|, Cˢ = 0.991 (κ = 0.41), half the surface TKE flux.\n"
     "GLORYS levels: 1 m top cell, 2 m at 10 m depth, 4 m at 25 m.",
     "scorecard_catke.png", mark_targets=True, separator=2.0)
print("wrote scorecard.png and scorecard_catke.png")
