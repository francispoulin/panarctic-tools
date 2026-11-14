# Patches

Numbered in the order to apply, across repositories. Each was checked to apply with `git am` on the head listed.
`TOOLS` below is this repository's clone. Before `git am`, check the branch with `git branch --show-current`.

| No. | Repository | Branch | Applies on | What it does |
|---|---|---|---|---|
| 001 | Oceananigans | `fjp/net-inflow-target` | 1e25007 | `NetInflow(T; depth_exponent)`: the net-inflow pin shifts each column ∝ Hᵖ; test of a common exterior datum |
| — | NumericalEarth | `nl+fjp/net-inflow-weighting` | — | already pushed (1a507d8): `net_inflow_depth_exponent` keyword |
| 002 | PanArcticShenanigans | `fjp/nested` | 3325594 | `[boundaries] net_inflow_depth_exponent`; NumericalEarth pin moved to 1a507d8; `runs/glorys-sqrtH-4d.toml` |
| 003 | PanArcticShenanigans | `fjp/nested` | after 002 | monitor: Flather's η term per edge, the pin's correction (Sv) and the largest velocity it adds |
| 004 | PanArcticShenanigans | `fjp/nested` | after 003 | GLORYS's mean η joined by a monotone cubic (no daily jumps); `runs/glorys-sqrtH-tendency-4d.toml` |
| 005 | PanArcticShenanigans | `fjp/nested` | after 004 | `gateways.jl`: hourly volume, heat and freshwater through Fram, BSO, Bering, Davis; `plot_gateways.jl` |
| 006 | PanArcticShenanigans | `fjp/nested` | after 005 | `runs/glorys-velocity-4d.toml`: the third weighting (p = 1) |
| 007 | PanArcticShenanigans | `fjp/nested` | ed899d3 | `gateways.csv` header: one name per column (fixes `plot_gateways.jl`) |
| 008 | PanArcticShenanigans | `fjp/nested` | d67e37d | `plot_gateways.jl` draws the hourly series as well as daily means |
| 009 | PanArcticShenanigans | `fjp/nested` | after 008 | `plot_gateways.jl`: reads the old header (no `sed`), whole days only, four-gateway sum as a fifth row |
| 010 | PanArcticShenanigans | `fjp/nested` | f8ff109 | `glorys_gateways`: GLORYS daily means through the same four sections; black dots on the gateway plot |
| 011 | PanArcticShenanigans | `fjp/nested` | after 010 | fixes 010 (GLORYS clock advanced before sampling); `setup.jl`; `glorys_gateways_run.jl` (CPU, no model) |

## Applying 001–006

    # Oceananigans clone
    git am $TOOLS/patches/001-*.patch && git push

    # PanArcticShenanigans, repository root
    git am $TOOLS/patches/00[2-6]-*.patch

    # PanArcticShenanigans/nested: fetch the new Oceananigans commit and NumericalEarth 1a507d8, then commit the Manifest
    julia +1.13 --project=. -e 'using Pkg; Pkg.update(["Oceananigans", "NumericalEarth"])'

## Runs after 006 (4 days, 1°, coupled, pressure off; the three pin weightings)

| Run file | Pin ∝ | Output |
|---|---|---|
| `runs/glorys-fastest-4d.toml` (existing; rerun for gateways.csv) | H⁰ (uniform per length) | `output-glorys-fastest-4d` |
| `runs/glorys-sqrtH-4d.toml` | √H | `output-glorys-sqrtH-4d` |
| `runs/glorys-velocity-4d.toml` | H | `output-glorys-velocity-4d` |

Then `runs/glorys-sqrtH-tendency-4d.toml` (smooth GLORYS tendency, 5-day relaxation).

Plots, from `nested/`: `julia --project=.. plot_budget.jl <folder>` and `julia --project=.. plot_gateways.jl <folder>`.

## What to look for

- `gateways.png`: daily means at Fram, BSO, Bering and Davis against the observed annual-mean ranges. The bands are
  annual means over many years and the run is four days of one January, so look for sign and order of magnitude.
- `nested_monitor.csv`: `pin_correction_Sv`, `pin_max_velocity` (how hard the pin pushes) and the Atlantic/Pacific
  split.
