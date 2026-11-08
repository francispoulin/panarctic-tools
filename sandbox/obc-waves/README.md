# obc-waves: which motions does an open boundary condition let out, and which does it keep?

A validation ladder for conditions on the barotropic transport V at an open edge, of the form
V − Vᵉˣᵗ = −α (ηᵇ − ηᵉˣᵗ) (Flather: α = c = √(gH); Dirichlet: α = 0). Linear rotating shallow water, f-plane,
open to the south, wall to the north. Fields vary as exp(ikx) along the edge, so each test is one-dimensional in y;
k = 0 is the edge-mean mode, the one that carries net volume.

| Rung | Test | Truth |
|---|---|---|
| 0 | steady balanced along-edge jet; external η correct, then 2 cm too high (a datum) | the jet, unchanged |
| 1 | one southward Poincaré packet, ω/f = 1.3, 1.6, 2, 3, 5 | r = (cₚ − α)/(cₚ + α), cₚ = c ω/√(ω² − f²) |
| 2 | a broadband pulse | the same r(ω) at every frequency |
| 3 | packets at 0°, 30°, 60° to the edge (k ≠ 0) | r = (Yₛ − α)/(Yₙ + α), `oblique.py` |

Results (`obc_1d_validation.png`): every run matches theory to within about 0.01, apart from a grid floor of
about lΔ/4 (ηᵇ is taken half a cell inside the edge).

- A wave's V/η is its phase speed cₚ, which grows without bound as ω → f. Flather's single c fits only ω ≫ f: it
  reflects 22 % at 1.3 f, 12 % at 1.6 f, 7 % at 2 f. A weight tuned to one frequency fixes that frequency only.
- At an angle Flather reflects far more: 29 % at 30° and 52 % at 60° for ω = 1.3 f.
- With a 2 cm datum, Flather builds a spurious boundary current of exactly gHδ/f (0.72 Sv here); in a closed basin
  that current runs along the walls and fills it (the obc-vortex E2b result).

## Rotating Flather (new, M7 in obc-vortex)

To second order in f²/ω², cₚ ≈ c/(1 − f²/2ω²). For k = 0 waves ∂u/∂t = f v, so this becomes a condition local in time:

    V = Vᵉˣᵗ − c (ηᵇ − ηᵉˣᵗ) + S/2,      dS/dt = −f H (uᵇ − uᵉˣᵗ),   S(0) = 0

- Reflection at normal incidence: 5 % at 1.3 f (Flather 22 %), 2 % at 1.6 f, about 1 % above 2 f.
- At 30°: 16 % at 1.3 f (Flather 29 %), 11 % at 2 f (14 %).
- The datum: S integrates the error in the along-edge velocity, which is datum-free, until the spurious current is
  gone. With η 2 cm too high the jet is recovered (0.004 cm rms at day 5, 0 by day 30; S/2c → −1.94 cm).
  Updating S once per 20-minute model step, as the Oceananigans callback does, gives the same result.
- Untested: an error in uᵉˣᵗ itself (S would then drive the model towards it), the 2D basin with Kelvin waves
  (`../obc-vortex/run_rotating.jl`), and the net transport.

## Files

| File | Contents |
|---|---|
| `sw1d.py` | the 1D solver (C-grid, RK4), theory, southward packets by spectral projection |
| `narrowband.py`, `broadband.py`, `steady.py` | rungs 1, 2, 0 |
| `oblique.py`, `oblique_rot.py` | rung 3, for Flather and rotating Flather |
| `rotating_flather.py` | rungs 0–1 for rotating Flather |
| `analytic_check.py` | checks the analytic packet used in `channel_waves.jl` |
| `figure.py` | `obc_1d_validation.png` |
| `channel_waves.jl`, `channel_plot.jl` | rungs 0–1 in Oceananigans (x-periodic channel), to check GravityWaveRadiation |

```bash
python3 narrowband.py; python3 broadband.py; python3 steady.py; python3 oblique.py
python3 rotating_flather.py; python3 oblique_rot.py; python3 figure.py         # about 30 minutes in all
julia -t 8 --project=$OCEAN_HPC_ROOT channel_waves.jl
julia --project=$OCEAN_HPC_ROOT/viz channel_plot.jl
```
