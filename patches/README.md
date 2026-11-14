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
| 012 | PanArcticShenanigans | `fjp/nested` | after 011 | `glorys_gateways_run.jl` always on the CPU (`PANARCTIC_ARCHITECTURE` overrides the config) |
| 013 | PanArcticShenanigans | `fjp/nested` | after 012 | `setup.jl` and `gateways.jl` import what they use (fixes `on_architecture not defined`) |
| 014 | PanArcticShenanigans | `fjp/nested` | 90cb4fd | `glorys_gateways_native.jl`: GLORYS gateway transports on GLORYS's own 1/12° grid (replaces `glorys_gateways_run.jl`); also through the model's staircase; `plot_gateways.jl` draws both |
| 015 | PanArcticShenanigans | `fjp/nested` | after 014 | fixes 014 (syntax error in `glorys_gateways_native.jl`) |
| 016 | PanArcticShenanigans | `fjp/nested` | after 015 | `runs/glorys-sqrtH-30d.toml` (1°, 30 days) and `runs/glorys-sqrtH-15km-5d.toml` (15 km, dt 10 min, ETOPO) |
| 017 | PanArcticShenanigans | `fjp/nested` | a342cff | GLORYS Arctic volume closure (`glorys_arctic_budget.csv`); `PANARCTIC_DOWNLOAD_ONLY=1`; `runs/tsubouchi-2005-15km-2mo.toml` (Sep–Oct 2005, 15 km) |
| 018 | PanArcticShenanigans | `fjp/nested` | 3ad3d86 | `plot_gateways.jl`: Tsubouchi et al. (2018) monthly volumes (orange), GLORYS Arctic budget panel, runs without a model; `nested/references/` |
| 019 | PanArcticShenanigans | `fjp/nested` | after 018 | fixes 017: longitude index in the full-longitude Arctic box (`point … outside the downloaded box`) |
| 020 | PanArcticShenanigans | `fjp/nested` | a7f3b53 | `[coupling] winds = false` (JRA55 winds zeroed); `runs/glorys-sqrtH-nowind-4d.toml` (1°, 4 days) |

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
