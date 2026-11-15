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

4. `julia --project=<PanArcticShenanigans>/nested ridge_fixes.jl [profile.csv] [days] [Nz]` — the ridge test with
   experimental fixes switched on and off in one session: partial cells (none, F1, F2, F1+F2) and whole cells (none,
   F2), f = 0 and f(72°N). Adds the background potential energy (water re-sorted by buoyancy into the model's own
   cells), so that ΔPE splits into mixing (ΔBPE) and conversion (ΔAPE = ΔPE − ΔBPE). Ends with a summary table.
   CPU; "none" should repeat `flow_over_ridge.jl`.

   The fixes are in `partial_cell_fixes.jl`, as method overrides loaded at run time (Oceananigans itself is not
   changed; see the file's header):
   - F1: in the z-star chain-rule term (∂x z)(∂z p′), take ∂z p′ with the full-cell face spacing, as p′ is
     integrated, instead of the partial spacing (partial cells only).
   - F2: immersed cells take the bottom wet cell's buoyancy in the p′ integration instead of T = S = 0 (both bottoms).
   With F1 on the run prints how many times the override was called; 0 means it never took effect.

Results (15 Nov, flow_over_ridge.jl): at 20 levels partial cells gain energy that whole cells don't (day 10, f = 0:
ΔE/KE₀ +0.26 vs +0.12; f(72°N): +0.09 vs +0.015 from day 2 to 10); at 50 levels partial ≈ whole.

**15 Nov, after `ridge_fixes.jl` and `probe_fixes.jl`:** the z-star chain-rule fixes F1/F2 work (F1 is called, it
changes the cliff force by 9e-8 m/s²) but change the 10-day energies by ≤ 0.001 KE₀. The pressure gradient is not
why partial cells do worse. The difference is all in the background potential energy: partial cells **mix** more
(f = 0: ΔBPE +0.41 vs +0.30 KE₀; f(72°N): +0.71 vs +0.48), with the fastest flow in the downstream bottom cell.

5. `julia --project=<PanArcticShenanigans>/nested mwe_partial_cells.jl [days] [Nz] [zstar|static]` — the minimal
   example for that: 2-D, f = 0, buoyancy tracer with b = Δb exp(z/d), whole vs partial cells, ΔBPE and a PASS/FAIL
   line (partial must not mix more than whole). This is the candidate test case for a PR.
