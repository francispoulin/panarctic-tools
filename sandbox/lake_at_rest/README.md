# Lake at rest over a cliff

Tests whether Oceananigans' horizontal pressure gradient is exactly zero for a horizontally uniform T(z), S(z) over a
steep step, as the 2005 pan-Arctic run had 3–5.6 m/s jets in bottom partial cells next to 650–1000 m steps.

1. `julia --project=<PanArcticShenanigans>/nested glorys_profile.jl` — writes `glorys_profile_0E_72N.csv` (GLORYS,
   0°E 72°N, 1 Sep 2005).
2. `julia --project=<PanArcticShenanigans>/nested lake_at_rest.jl [profile.csv] [Nz]` — CPU, a few minutes.

Part 1 (t = 0): the force −∂x p′ at every active u face, split into the part along the levels and the z-star chain-rule
part, for partial vs whole cells and η = 0 vs −0.47 m; then again with the immersed cells given the column's bottom
T and S instead of 0. Part 2: two days of stepping, largest |u| every 6 h.

3. `julia --project=<PanArcticShenanigans>/nested flow_over_ridge.jl [profile.csv] [days]` — stratified flow over a
   periodic ridge (two 1000 m cliffs) with an energy budget (kinetic, potential, free surface; linear equation of
   state), partial vs whole cells × 20 vs 50 levels × f = 0 vs 72°N, ten days each. Energy should only fall; a rise,
   or fast flow in bottom cells, points at the partial-cell bottom. CPU, some minutes per case.

Results so far (14 Nov): the lake at rest is exact at η = 0 for partial and whole cells, 20 and 50 levels; at
η = −0.47 m both have a ~1.5–2 × 10⁻⁷ m/s² z-star datum error (u stays ~2–3 mm/s). Filling the immersed cells changes
nothing. So the static pressure gradient is not the cause of the deep jets.
