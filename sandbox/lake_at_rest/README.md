# Lake at rest over a cliff

Tests whether Oceananigans' horizontal pressure gradient is exactly zero for a horizontally uniform T(z), S(z) over a
steep step, as the 2005 pan-Arctic run had 3–5.6 m/s jets in bottom partial cells next to 650–1000 m steps.

1. `julia --project=<PanArcticShenanigans>/nested glorys_profile.jl` — writes `glorys_profile_0E_72N.csv` (GLORYS,
   0°E 72°N, 1 Sep 2005).
2. `julia --project=<PanArcticShenanigans>/nested lake_at_rest.jl [profile.csv] [Nz]` — CPU, a few minutes.

Part 1 (t = 0): the force −∂x p′ at every active u face, split into the part along the levels and the z-star chain-rule
part, for partial vs whole cells and η = 0 vs −0.47 m; then again with the immersed cells given the column's bottom
T and S instead of 0. Part 2: two days of stepping, largest |u| every 6 h.
