# CATKE under sea ice: single columns against LES

Plan: project doc `claude/catke-les-plan.md`. Truth: 3-D large-eddy simulations (Oceananigans `NonhydrostaticModel`,
~1 m, implicit LES with WENO order 9, as LESbrary.jl and Oceananigans' Langmuir example); tested: hydrostatic single
columns with CATKE variants at 1, 4 and 8 m, with the same forcing and initial state.

- `cases.jl`: the cases. Step 0 (now): two cases of CATKE's own calibration (LESbrary's 24-hour suite): free
  convection, strong wind, and strong wind with Stokes drift (Langmuir). Stage 1 (sea ice, GLORYS September profile,
  TEOS-10, quadratic ice–ocean drag) follows once step 0 checks out.
- `les.jl <case> [small]`: the LES (GPU). Default 256 × 256 × 128 m at 1 m; `small` = 128 × 128 × 128 m at 2 m.
- `columns.jl <case>`: CATKE columns (CPU, a minute or two): stock, and without the shear-driven surface TKE flux.
- `compare.jl <case> [small]`: profiles at 6, 12, 24 h and a table of errors (`output/<case>_comparison.md`, `.png`).

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
