# Net-inflow pin weighting, smooth GLORYS target, pin monitor (2026-11-13)

Each patch was checked to apply with `git am` on the head named below.

| Repository | Branch, head it applies on | Patches |
|---|---|---|
| Oceananigans | `fjp/net-inflow-target`, 1e25007 | `oceananigans/0001` — `NetInflow(T; depth_exponent)` + test |
| NumericalEarth | `nl+fjp/net-inflow-weighting`, 1a507d8 | none: already pushed (`net_inflow_depth_exponent` keyword) |
| PanArcticShenanigans | `fjp/nested`, 3325594 | `panarcticshenanigans/0001–0003` |

## Steps

`TOOLS` is this folder's clone (`panarctic-tools`); run each block from the repository named.

1. Oceananigans (your clone with `fjp/net-inflow-target`):

       git branch --show-current          # must print fjp/net-inflow-target
       git am $TOOLS/patches/2026-11-13-pin-weighting/oceananigans/*.patch
       git push

2. PanArcticShenanigans (repository root):

       git branch --show-current          # must print fjp/nested
       git am $TOOLS/patches/2026-11-13-pin-weighting/panarcticshenanigans/*.patch

3. Environment, from `nested/` (fetches the new Oceananigans commit and NumericalEarth 1a507d8):

       julia +1.13 --project=. -e 'using Pkg; Pkg.update(["Oceananigans", "NumericalEarth"])'

   Then commit the changed `nested/Manifest.toml` and push.

4. Runs, in this order (each 4 days at 1°, coupled, pressure off):

   - `runs/glorys-sqrtH-4d.toml` — only change from `glorys-fastest-4d`: pin ∝ √H.
   - `runs/glorys-sqrtH-tendency-4d.toml` — adds the (now smooth) GLORYS tendency, relaxation 5 days.

   Plots as before: `julia --project=.. plot_budget.jl <folder>`.

## What to look for

- `pin_correction_Sv` and `pin_max_velocity` in `nested_monitor.csv` (also on the hourly log line): how hard the pin
  pushes. Both runs log it; compare with `output-glorys-fastest-4d` re-run if you want the p = 0 numbers.
- `atlantic_in_Sv`, `pacific_in_Sv`: the split between the edges. The sea-level datum offset should no longer
  drive a shelf circulation along the edges with √H; the Pacific mean inflow (0.24 Sv with p = 0) is the number to
  watch.
- In the tendency run: net inflow without daily jumps, and mean η on GLORYS's with no one-day lag.
