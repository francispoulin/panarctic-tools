# CATKE under sea ice: single columns against LES

Plan: project doc `claude/catke-les-plan.md`. Truth: 3-D large-eddy simulations (Oceananigans `NonhydrostaticModel`,
~1 m, implicit LES with WENO order 9, as LESbrary.jl and Oceananigans' Langmuir example); tested: hydrostatic single
columns with CATKE variants, k-ε and NEMO's TKE at 1, 4 and 8 m, with the same forcing and initial state.

- `cases.jl`: the cases: LESbrary's 24-hour suite (free convection, weak wind strong cooling, medium wind medium
  cooling, strong wind weak cooling, strong wind), each wind case without waves and with the calibration's broadband
  Stokes drift (`<case>_stokes`, LESbrary's own case). Stage 1 (sea ice) follows. Also `column_profiles(file)`, which
  reads both kinds of column output.
- `les.jl <case> [small|fine]`: the LES (GPU), 256 m deep. Default 512 × 512 m at 2 m (the calibration extent);
  `small` = 256 × 256 m at 2 m; `fine` = 256 × 256 m at 1 m.
- `columns.jl <case>`: CATKE columns (CPU, a minute or two): stock; without the shear-driven surface TKE flux
  (Cᵂu★ = 0); rigid wall (Cᵂu★ = 0 and Cˢ = 0.967, von Kármán constant 0.4 instead of the wave-influenced 0.47); k-ε.
- `nemo_tke.jl <case>`: NEMO's TKE scheme (zdftke.F90, reference namelist, with enhanced vertical diffusion for
  convection) as a standalone column (CPU, seconds): `nemo` (Langmuir cells and TKE penetration on), `nemo_no_waves`,
  `nemo_lc_only`, at 1, 4, 8 m. Picked up by compare.jl and compare_lesbrary.jl (dash-dot lines).
- `run_suite.jl [cases...]`: LES (skipped if done), columns, NEMO columns and comparisons for several cases.
  `run_suite.jl rung1`: rung 1 of the validation ladder (plan doc): five free-convection cases into uniform N²
  (entrainment law h² N² = 2(1 + 2A) Jb t) and `arctic_brine_convection` (brine rejection at freeze-up, TEOS-10,
  weak mixed layer over a strong halocline, 4 days), then `scaling.jl` (tables and pass/fail per column).
  `run_suite.jl rung3`: stress without waves: Kato–Phillips (Price's law), rotating stress, and drifting ice
  (quadratic drag) over the Arctic profile alone, with brine and with melt; then `stress_scaling.jl`.
- `calm_melt.jl`: rung 2, melt under calm conditions (columns only, CPU): does any scheme mix the freshwater
  below the top cell when nothing drives turbulence?
- Columns run at 1, 4, 8 m and on the pan-Arctic model's vertical grids (`nz50`, `nz100`: setup.jl's exponential
  grid cut at the domain depth).
- `scan.jl <case> [Δz]`: step 2, a scan of Cˢ × Cᵂu★ in columns against our LES (default strong_wind, no waves).
- `compare_lesbrary.jl <their case> [hours]`: step 1, our LES and columns against the calibration LES (LESbrary
  profiles at 1, 2, 4 m from glwagner/SingleColumnModelCalibration.jl, downloaded to `data/profiles/`; see the README
  there for the Dropbox link). Their `strong_wind` has Stokes drift: compared with our `strong_wind_stokes`.
- `compare.jl <case> [small|fine]`: profiles at 6, 12, 24 h and a table of errors (`output/<case>_comparison.md`, `.png`).

Run in a REPL from the PanArcticShenanigans clone's environment (same Oceananigans as the model):

```
cd /mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-061
julia --project=.
```
```julia
using CUDA; CUDA.device!(3)
empty!(ARGS); push!(ARGS, "free_convection", "small")
include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/les.jl")
empty!(ARGS); push!(ARGS, "free_convection")
include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/columns.jl")
empty!(ARGS); push!(ARGS, "free_convection", "small")
include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/compare.jl")
```

Output goes to `les_under_ice/output/` (not committed).
