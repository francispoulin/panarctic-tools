# obc-vortex: open boundary conditions for the barotropic transport

A minimal test of open boundary conditions for the depth-integrated (barotropic) transport, built to understand why
the pan-Arctic model lets 10–15 Sv too much water in through its Atlantic edge (October–November 2026).

Linear rotating shallow water: Oceananigans `HydrostaticFreeSurfaceModel` with Nz = 1, `SplitExplicitFreeSurface`,
an f-plane and no advection, so the open-boundary machinery is the same as in the pan-Arctic model. A 2000 km square
basin with walls on the west, north and east and an open south side. A geostrophic Gaussian vortex centred on the open
side crosses it inwards on one half and outwards on the other, with zero net transport: an exact steady solution, so
the right answer is known. Theory, notation and the derivation of each condition: `balanced_obc_notes.tex`.

Written and run with Oceananigans 0.113.6 on CPU (about 3 minutes for E0–E2 with every method); updated for
Oceananigans 0.117 (oceananigans-hpc v1.4): the open south side is declared `Open` in the grid topology
(`TwoSided{Open, Closed}` in y), as open-boundary schemes now require.

On oceananigans-hpc, run the simulations in the platform environment and the plots in `viz/`:

```bash
julia -t 8 --project=$OCEAN_HPC_ROOT run_davies.jl
julia --project=$OCEAN_HPC_ROOT/viz compare.jl davies
```

The older drivers (`run_all.jl`, `run_m6.jl`, `run_e2b.jl`, `run_trust.jl`) run and plot in one go, so they need an
environment with both Oceananigans and CairoMakie.

## Files

| File | Contents |
|---|---|
| `obc_vortex.jl` | parameters, grid, vortex, the boundary conditions, `run_case`, `run_reference` |
| `pv_inversion.jl` | PV inversion on the C-grid (sparse LU), balanced/unbalanced split, energy |
| `analysis.jl` | `compare(experiment)`: net transport, unbalanced energy, error against the truth, mean η; maps |
| `run_all.jl` | every method on E0, E2, E1 (or the experiments given), the reference where needed, the comparisons |
| `run_m6.jl` | M6 with T_C = 12 h, 2 d, 5 d on E0, E2, E1 |
| `run_e2b.jl` | the growing offset (E2b) |
| `run_trust.jl` | `shelf`: E1, E2, E2b with a shelf (E4); `net_inflow`: E3, 30 days |
| `run_davies.jl` | M3 (Dirichlet + sponge) with rims of 200, 500, 1000 km at 2 h and 1000 km at 6 h, and M1 and M5b each with a 200 or 500 km sponge, on E1 and E2b |
| `compare.jl` | figures and tables from `results/` in a plotting-only environment (`compare.jl davies`, or experiment names) |
| `../balanced-flather-block.jl` | M5b written for the pan-Arctic model (Atlantic and Bering edges) |

```bash
cd sandbox/obc-vortex
julia --project=<env with Oceananigans and CairoMakie> run_all.jl          # results/*.jls, figures/*.png
julia --project=... run_m6.jl
julia --project=... run_e2b.jl
julia --project=... run_trust.jl shelf
julia --project=... run_trust.jl net_inflow
```

## Methods (condition on the barotropic transport V at the open faces)

| | |
|---|---|
| M1 `:flather` | V = Vᵉˣᵗ − c (ηᵇ − ηᵉˣᵗ), c = √(gH) |
| M2 `:dirichlet` | V = Vᵉˣᵗ |
| M3 `:davies` | M2 plus relaxation of u, v to the external state in a 200 km rim |
| M4 `:flather_target` | M1 with `target_transport` (constant; uniform shift per unit length) |
| M5b `:balanced_net` | V = Vᵉˣᵗ − c (ηᵇ − ηᵉˣᵗ − C), C chosen so the net transport equals the external net (√(gH)-weighted) |
| M5a `:balanced_pv` | same, C from the PV inversion with a volume constraint |
| M6 `:balanced_filtered` | same, C a slow average (T_C) of the √(gH)-weighted mean edge mismatch |

## Experiments

E0 `:steady`; E1 `:adjustment` (unbalanced bump, waves must leave; reference: basin extended 8000 km south);
E2 `:offset` (external η 3 cm higher from t = 0); E2b `:offset_ramp` (offset grows over 1 day);
E3 `:net_inflow` (1 Sv, 10-day period, 30 days); E4 `shelf = true` (50 m shelf along the west of the open side).

Defaults: L = 2000 km, Δ = 20 km, H = 500 m (c = 70 m/s, L_d = 500 km), f = 1.4e-4 s⁻¹, vortex R = 150 km,
η₀ = 0.1 m; bump 0.2 m, 150 km at (1000, 1500) km.

## Results (rms η error against the truth, cm)

| | M1 | M2 | M3 | M4 | M5b | M5a | M6 (2 d) |
|---|---|---|---|---|---|---|---|
| E1 waves | **0.05** | 0.46 | 0.38 | 0.32 | 0.32 | 0.25 | **0.05** |
| E2 offset at t = 0 | 2.13 | **0** | **0** | **0** | **0** | 0.03 | **0** (initialised with it) |
| E2b growing offset | 2.14 | **0** | **0** | — | **0** | 0.07 | 1.92 |
| E1 shelf | 0.07 | 0.42 | 0.35 | 0.29 | 0.28 | 0.25 | — |
| E2b shelf | 2.34 | **0** | **0** | 0.11, shelf jets | **0** | 0.07 | — |
| E3 net inflow | 35 % too little | exact | exact | none (constant target) | exact | 60 % too much, drifts | — |

- Flather does what it is designed to: its transport error is c × (mean edge η mismatch) × width. Best for letting
  waves out, wrong whenever the external and model free surfaces differ by a datum.
- There is no time-scale gap for the barotropic mode: the basin fills to the external η in about L/c, the same time
  waves take to leave, so a time-filtered C (M6) only works if the offset is known at t = 0.
- M5a assumes the balanced part keeps all the basin's volume; waves carry volume out of an open basin, and under a
  net inflow the assumption feeds back on itself. Dropped.
- Decision: trust the transport at low frequency. M5b is the best transport-trusting method: exact net transport,
  immune to datum offsets, no shelf jets, and fewer trapped waves than Dirichlet. Its cost: the part of an outgoing
  wave that is uniform along the edge (the net-volume part) reflects.

## Notes

- M2/M3 use a small callable (`DirichletTransport`) because a plain normal-flow condition is also filled without a
  clock at initialisation, which `DiscreteBoundaryFunction` does not accept.
- M5b computes C inside the boundary function at every barotropic substep; in Oceananigans it would be a
  time-dependent, √(gH)-weighted `target_transport`.
- The energy split has a floor of about 6e7 J m³/kg from discretisation.
- Oceananigans ≥ 0.117: the basin's y topology is `TwoSided{Open, Closed}`; the reference basin stays `Bounded`.
