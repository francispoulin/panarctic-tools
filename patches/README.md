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
| 021 | PanArcticShenanigans | `fjp/nested` | a049a53 | `jra55_read_fix.jl`: JRA55 reads for runs not starting on 1 January (fixes `-(::Nothing, ::Int64)` in `build_flux_routing`); NumericalEarth bug, to report upstream |
| 022 | PanArcticShenanigans | `fjp/nested` | after 021 | `[grid] bottom = "full"` (GridFittedBottom); run files `tsubouchi-2005-15km-fullcell-10d`, `tsubouchi-2005-15km-nz50-10d` |
| 023 | PanArcticShenanigans | `fjp/nested` | 7063a25 | `runs/tsubouchi-2005-15km-fullcell-2mo.toml` (Sep–Oct 2005, 15 km, whole cells) |
| 024 | PanArcticShenanigans | `fjp/nested` | after 023 | `runs/tsubouchi-2005-15km-fullcell-nz50-2mo.toml` (Sep–Oct 2005, 15 km, whole cells, 50 levels) |
| 025 | PanArcticShenanigans | `fjp/nested` | 3cdb26e | Nares, Jones, Barrow sections in `gateways.jl` (GLORYS budget still the four); daily ice h, ℵ to `nested_ice.jld2`; CAA summary in `plot_gateways.jl`; `runs/tsubouchi-2005-15km-fullcell-caa-2mo.toml` |
| 026 | PanArcticShenanigans | `fjp/nested` | after 025 | `runs/tsubouchi-2005-15km-partial-{hybrid,ke}-15d.toml`: partial cells, 15 days, with a momentum fix from `sandbox/lake_at_rest/thickness_jump_fixes.jl` loaded first (see the files' headers) |
| 027 | PanArcticShenanigans | `fjp/nested` | after 026 | `heat_freshwater_tsubouchi.jl` (gateway heat and freshwater on Tsubouchi et al.'s references, model and GLORYS vs their Table 3); `references/tsubouchi2018_table3.csv` |
| 028 | PanArcticShenanigans | `fjp/nested` | after 027 | `heat_freshwater_tsubouchi.jl` compares Tsubouchi et al.'s monthly blocks (same days) for volume, heat and freshwater; `references/tsubouchi2018_monthly_{heat,freshwater}.csv` (from Figs. 7a, 9a) |
| 029 | PanArcticShenanigans | `fjp/nested` | after 028 | fix applied by hand (no patch file): in `plot_gateways.jl` line 181, `"$a–$b"` → `"$(a)–$(b)"` (parse error from 025) |
| 030 | Oceananigans | new branch `fjp/mass-consistent-bottom-advection` from `main` (9840bec, v0.118.1) | — | mass-consistent momentum advection at bottom-cell thickness jumps (`src/Advection/mass_consistent_bottom_advection.jl`) and its test (`test/hydrostatic_free_surface/`); the fix of `sandbox/lake_at_rest/thickness_jump_fixes.jl` (hybrid, wet cells only) as source |
| 031 | PanArcticShenanigans | `fjp/nested` | after 028 | `[output] ocean_3d_days` (daily averages of 3-D T, S, u, v in `nested_ocean_3d.jld2`); `runs/tsubouchi-2005-15km-partial-hybridwet-2mo.toml` (production settings, 61 days, fix with wet cells only) |
| 032 | PanArcticShenanigans | `fjp/nested` | after 031 | `section_properties.jl` (model inflow/outflow and mean sections from the 3-D output), the same for GLORYS in `glorys_gateways_native.jl`, `section_properties_common.jl`, `plot_sections.jl` (heat and freshwater errors split into volume and property parts; sections model vs GLORYS) |
| 033 | Oceananigans | `fjp/mass-consistent-bottom-advection` (PR #6533) | after 030 | the fix only for vector-invariant schemes whose kinetic-energy-gradient scheme has momentum fluxes (WENO, centred, upwind); `VectorInvariant()` (EnergyConserving) keeps the stock form (CI failure in `simulation/netcdf_writer`) |
| 034 | PanArcticShenanigans | `fjp/nested` | after 032 | inflow/outflow and sections per 30 km stretch of the line and level (`binned` in `section_properties_common.jl`), not per staircase face; rerun `glorys_gateways_native.jl`, `section_properties.jl`, `plot_sections.jl` |
| 035 | PanArcticShenanigans | `fjp/nested` | after 034 | `rivers_and_straits.jl` (CPU): does the Baltic connect to the Skagerrak, map of the Danish straits, JRA55 river discharge by region against what the routing delivers |
| 036 | PanArcticShenanigans | `fjp/nested` | after 035 | `rivers_and_straits.jl` loads `jra55_read_fix.jl` (the river part failed in `jra55_read_data`) |
| 037 | PanArcticShenanigans | `fjp/nested` | after 036 | `[grid] bathymetry = "GLORYS"` (`glorys_bathymetry.jl`, GLORYS's `deptho` regridded like `regrid_bathymetry`); `runs/tsubouchi-2005-15km-partial-hybridwet-nz50-glorysbathy-2mo.toml` |
| 038 | PanArcticShenanigans | `fjp/nested` | after 037 | `coastal_salinity.jl`: top-20 m salinity at seven stations from the Baltic to the BSO coast, model (3-D output) against GLORYS (parent files), day by day |
| 039 | PanArcticShenanigans | `fjp/nested` | after 038 | `plot_sections.jl` and `coastal_salinity.jl` average GLORYS over each model averaging window (needed for 2-day averages, e.g. the 10 km run) |
| 040 | PanArcticShenanigans | `fjp/nested` | after 039 | `coastal_salinity.jl` also gives the freshwater content of the top 200 m at each station (m of fresh water relative to 34.8), model and GLORYS |
| 041 | Oceananigans | `fjp/net-inflow-target` (PR #6524) | after 53efe85 | net-inflow pin: no local variable across `@synchronize` (CPU `UndefVarError: Ny`); GPU unaffected |
| 042 | PanArcticShenanigans | `fjp/merge-nested` (new, from `main` 14427af) | first | `config_smoke_test.toml` (~100 km, 1–5 Jan 1993) and `analysis/compare_transports.jl` (largest difference of every transport between two experiments) |
| 043 | PanArcticShenanigans | `fjp/merge-nested` | after 042 | pins NumericalEarth → `fjp/net-inflow` e22c1bfd (#938 on #676 head b2eb617d) and Oceananigans → `fjp/net-inflow-target` 705b9e7c (#6524); defaults unchanged, so the smoke test should match 042 |
| 044 | PanArcticShenanigans | `fjp/fix-online-transports` (new, from `main` 14427af; own PR to main) | first | online transports after #43: child tracers, child grid, writer on the child (was `ocean.model.tracers`: FieldError, every run stopped); `$(λ)°E` in the section-end warning. `fjp/merge-nested` is rebased on top of it |
| 045 | PanArcticShenanigans | `fjp/merge-nested` | after 043 (on top of 044) | volume control on by default (`volume_control.jl`, from `glorys_volume_target.jl`): net inflow follows GLORYS's domain-mean η, relaxation 1 day, pin ∝ H^½; `net_inflow_Sv` or Flather alone as alternatives; hourly volume log line |
| 046 | PanArcticShenanigans | `fjp/merge-nested` | after 045 | hourly `monitor.csv` (`[output] monitor_hours = 1`): open-side transports, mean η vs GLORYS, volume budget, fastest point; `plotting/volume_budget.jl`, `plotting/max_speed.jl` (from fjp/nested, without the dropped volume-control diagnostics) |
| 047 | PanArcticShenanigans | `fjp/merge-nested` | after 046 | `[atmosphere] inverse_barometer = true` (default; #676 option, issue #37); `jra55_read_fix.jl` (runs can start on any date); `config_smoke_test_sep2005.toml` |
| 048 | PanArcticShenanigans | `fjp/merge-nested` | after 047 | `partial_cell_fix.jl` (Oceananigans #6533 as run-time override, `[grid] partial_cell_fix = true`); `[grid] bathymetry` GEBCO2026/ETOPO/GLORYS (`glorys_bathymetry.jl`), `[grid] bottom` partial/full |
| 049 | PanArcticShenanigans | `fjp/merge-nested` | after 048 | `nested_pressure_fix.jl`: `forcing_barotropic_potential(::NestedModel)` → child; without it no nested run (main since #43) gets the JRA55 pressure, and the inverse barometer is zero |
| 050 | PanArcticShenanigans | `fjp/test-050` only (diagnostic, not for the merge) | after 049 | `config_smoke_test_full_cells.toml`: smoke test with whole bottom cells, for the fast v at (25, 1, 4) in the partial bottom cell of the Atlantic open side |
| 051 | PanArcticShenanigans | `fjp/merge-nested` | after 049 | `analysis/glorys_transports.jl` (GLORYS's daily transports through the [transports] sections on GLORYS's own grid, with transports.jl's own code), `plotting/gateway_transports.jl` (model vs GLORYS) |
| 052 | PanArcticShenanigans | `fjp/test-052` only (diagnostic, not for the merge) | after 051 | `config_cost_fix_on.toml`, `config_cost_fix_off.toml`: 15 km, 1 day, partial-cell fix on/off, monitor and online transports off, for the cost question on Oceananigans #6533 |
| 053 | PanArcticShenanigans | `fjp/test-053` only (diagnostic, not for the merge) | after 051 | `boundary_wet_center_fix.jl` (candidate #676 fix: GLORYS sampled at the centre of the wet part of partial bottom cells; `[boundaries] sample_wet_center = true`), `config_smoke_test_wet_center.toml`, `analysis/boundary_bottom_cells.jl` (GLORYS at the whole-cell vs wet centre on the boundary rows, with the model's value) |
| 054 | PanArcticShenanigans | `fjp/test-054` only (diagnostic, not for the merge) | after 051 | `[atmosphere] pressure = false` skips `nested_pressure_fix.jl`; `config_pressure_on.toml`, `config_pressure_off.toml` (15 km, 4 days from 1 Jan 1993, pressure + inverse barometer both on / both off) to compare with GLORYS |
| 055 | PanArcticShenanigans | `fjp/test-053` only (diagnostic) | after 053 | `analysis/boundary_bottom_cells.jl`: KernelAbstractions through Oceananigans (it is not a dependency of the project) |
| 056 | PanArcticShenanigans | `fjp/test-053` only (diagnostic) | after 055 | the Atlantic face jet at (25, 1, 4): `config_smoke_test_inflow_1h.toml` (inflow nudging 1 h) and `config_smoke_test_face_direction.toml` + `radiation_face_direction.jl` (radiation decides inflow/outflow on the boundary face, not the first interior face) |

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
