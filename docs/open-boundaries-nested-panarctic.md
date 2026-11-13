# Open boundary conditions of the nested pan-Arctic model

*Francis Poulin, 13 November 2026. The model: PanArcticShenanigans `nested/`, built on NumericalEarth's
`nested_ocean_model` (PR #676) with our additions (PR #923, branch `nl+fjp/nested-panarctic`; `net_inflow` on branch
`nl+fjp/net-inflow`), Oceananigans 0.117.2 plus the `NetInflow` target (branch `fjp/net-inflow-target`).*

This note describes what the model imposes at its two open edges, with equations, and how the default set-up fails
to control the domain's volume. It then describes the fix we now use (sections 4.4 and 4.5) and how we are testing it.

## 1. Set-up and notation

- **Domain.** Conformal conic grid centered on the pole (λ₀ = −40°), 1° (100 km) for tests and 15 km for production, 20
  levels, z-star vertical coordinate. Two open edges: the **Atlantic** edge (south side of the grid, faces j = 1) and
  the **Pacific/Bering** edge (north side, faces j = Ny + 1). The other sides are land or closed.
- **Parent.** GLORYS12 daily means, interpolated in time between consecutive days and in space onto the edge cells.
  Velocities are rotated into the grid's frame.
- **Free surface.** Split-explicit. The depth-integrated (barotropic) flow and η are advanced in about 30 short
  substeps within each model time step; the depth-dependent (baroclinic) flow, T and S are advanced once per step.

Notation at one open edge (normal velocity v; the same holds for u on an east or west edge):

| Symbol | Meaning |
|---|---|
| $v^G(z)$, $\eta^G$, $T^G$, $S^G$ | GLORYS at the edge, interpolated in time |
| $\bar v^G$ | GLORYS's depth mean of $v^G$ over GLORYS's own water column |
| $H$ | the model's wet depth at the edge face |
| $c = \sqrt{gH}$ | barotropic gravity-wave speed |
| $V$ | the model's barotropic (depth-integrated) transport per unit width, m² s⁻¹ |
| $\eta_b$ | the model's η at the edge face (average of the edge cell and its halo) |
| $s$ | +1 on the south (Atlantic) edge, −1 on the north (Pacific) edge, so that $s\,V$ is the inflow |
| $\Delta x$ | width of a face along the edge |
| $Q_i = s_i \sum V \Delta x$ | net inflow through edge $i$, m³ s⁻¹ (Sv = 10⁶ m³ s⁻¹) |

## 2. What #676 imposes at each open edge

Each open edge gets five conditions, plus a relaxation band and a bathymetry blend just inside it.

### 2.1 Barotropic transport: Flather, at every substep

The exterior transport is GLORYS's depth mean times the *model's* wet depth:

$$V^{ext} = \bar v^G H .$$

At every barotropic substep, after the momentum update, the transport at the edge face is overwritten by Flather's
condition:

$$V_b = V^{ext} - s\, c \,(\eta_b - \eta^G).$$

When the model's sea level at the edge stands above GLORYS's, water flows out; when it stands below, water flows in.
**This is the flux that changes η and hence the volume.** It is the only route by which GLORYS's η reaches the model.

### 2.2 Free surface at the edge: Chapman radiation

The edge value of η radiates outward at speed $c$ (implicit form, as in ROMS):

$$\eta_b^{n+1} = \frac{\eta_b^{n} + C\,\eta_1^{n+1}}{1 + C}, \qquad C = \frac{c\,\Delta t}{\Delta y},$$

where $\eta_1$ is the first interior cell. It sets η in the halo just outside the edge, which enters only two places:
the pressure gradient on the edge face in the barotropic momentum equation, and Flather's $\eta_b$ (the average of the
edge cell and this halo value). The η of the interior cells changes only through the divergence of the transport. It
does not relax η toward GLORYS.

### 2.3 Full-depth normal velocity: "consistent" profile with oblique radiation

The external profile is GLORYS's, shifted uniformly so that it integrates to $V^{ext}$ over the *model's* column (the
shift is zero where the model's column at the edge matches GLORYS's):

$$v^{ext}(z) = v^G(z) + \frac{\bar v^G H - \sum_k v^G_k \Delta z_k}{H}, \qquad \sum_k v^{ext}_k \Delta z_k = \bar v^G H = V^{ext}.$$

The model's velocity at the edge face follows an oblique radiation scheme: on inflow it is nudged toward
$v^{ext}(z)$ with a 9-hour time scale; on outflow it radiates from the interior. **These velocities carry T and S
across the edge.**

### 2.4 Tangential velocity and tracers

- Tangential velocity: the same oblique radiation scheme as the normal velocity, with GLORYS's tangential velocity as
  the external value (nudged toward it over 9 hours on inflow, radiated on outflow). No depth-mean adjustment.
- T and S: a tracer reservoir (20 km): inflowing water carries $T^G$, $S^G$; outflowing water carries the model's
  values.

### 2.5 Relaxation band (sponge)

Within 2/3° of each open edge, T, S and the baroclinic part of the velocity (the profile minus its depth mean) are
relaxed toward GLORYS:

$$\partial_t \psi = \dots - \mu(x)\,(\psi - \psi^G), \qquad \mu = 1/(6\,\text{h}) \text{ at the edge, decreasing inward.}$$

The depth-mean flow and η are **not** relaxed.

### 2.6 Bathymetry blend

Within 1/3° of each open edge the model's seafloor is blended toward GLORYS's, so that the columns at the edge match
those GLORYS's state was computed over.

### 2.7 Initial state

GLORYS interpolated onto the grid at the start date. It is not dynamically balanced.

## 3. How the default set-up fails

### 3.1 Flather sets the volume through the sea-level mismatch, not through GLORYS's transport

Summing Flather over an edge gives the net inflow

$$Q = \sum s\, \bar v^G H\, \Delta x \;+\; \sum c\,(\eta^G - \eta_b)\,\Delta x .$$

The second term is very stiff. For our edges, $\sum c\,\Delta x$ is about 8 × 10⁸ m² s⁻¹: **a mismatch of 1 cm between the
model's and GLORYS's sea level along the edges drives about 8 Sv.** With a domain area $A \approx 2.1\times10^{13}$ m²,
the domain-mean η relaxes to the edge value of GLORYS on a time scale

$$\tau = \frac{A}{\sum c\,\Delta x} \approx 7 \text{ hours}.$$

So the domain's volume is set by whatever keeps the model's edge sea level at GLORYS's. When wind changes the model's
sea level inside the domain relative to its edges differently from GLORYS, the whole interior moves up or down with it.
In a 30-day coupled run (1°, January 1993, pressure off) this produced swings in domain-mean η of +11 cm and −10 cm
relative to GLORYS, all of it through the edges (volume budget, section 5).

### 3.2 GLORYS's transport, interpolated onto our edges, does not close

$V^{ext} = \bar v^G H$ uses the model's depth, not GLORYS's, and GLORYS's velocities interpolated onto our faces. At
15 km, the Atlantic edge's net $\sum \bar v^G H\,\Delta x$ is −6.8 Sv against −0.7 Sv computed on GLORYS's own grid.
The gross flow is about 85 Sv each way, so a 4% error in sampling or depth gives several Sv of net. Even on its own
grid, GLORYS's edge transports do not match its sea-level change: over January 1993 they imply a 35 cm fall in mean
η, while GLORYS's mean η fell 3 cm. GLORYS assimilates sea level, and its transports are not consistent with it.

### 3.3 The water and the velocities carrying T and S (checked: they agree)

From 2.1 and 2.3 alone one would expect the full-depth velocities at the edge faces to integrate to $\bar v^G H$ and the
water to cross with $V_b$, a difference of Flather's η term. The 30-day run with plain #676 shows otherwise: the net
inflow carried by the velocities at the edge faces matches the inflow that changes the volume to within a few tenths
of a Sv all month (section 5). The model's barotropic correction brings the edge velocities into line with $V_b$.
This is not a failure of the default set-up.

### 3.4 The stored averaged transport is not the flux that moved the water

The averaged barotropic transport stored at the edge faces (used for tracers) is overwritten once at the end of each
step by Flather with the final η. It is not the average of the substep values that changed η. Without a constraint
the two are close; with an outside correction they diverged by tens of Sv (section 5).

## 4. What we added, and the fix

### 4.1 `zero_net_inflow` (PR #923, first commit)

At every new GLORYS sample, one depth-uniform velocity is removed from $\bar v^G$ on all open edges:

$$\bar v^G \leftarrow \bar v^G - s_i \frac{Q^{ext}}{A_b}, \qquad Q^{ext} = \sum_i s_i \sum \bar v^G H\,\Delta x, \quad A_b = \sum_i \sum H\,\Delta x .$$

This removes the net of the first term in 3.1, which fixes 3.2 for the total. It does nothing about the η term. In 30
days of coupled runs it reduced the error in mean η from −9.7 to −2.1 cm at day 30 (peak +11 → +4 cm), but the
swings from 3.1 remain.

### 4.2 Prototype controls from outside the substeps (not kept)

We tried holding the net inflow at zero from outside the model, once per time step, by shifting either $\bar v^G$ or
$\eta^G$ on both edges.

- **Measured with the stored face transport,** the control drove that number to zero while the volume still swung
  as before (3.4).
- **Measured with the volume,** the volume was held, but the shift of $\bar v^G$ needed swung between −140 and +100 Sv,
  and the velocities carrying T and S (2.3) then carried ±40 Sv of net flow that never entered the volume (3.3).

Correcting the data Flather reads, from outside its substeps, is the wrong lever.

### 4.3 The fix: pin each edge's inflow inside the substeps (PR #923, `edge_transports`)

This is the volume constraint of NEMO (`ln_vol`) and ROMS (`VolCons`), applied edge by edge. Given a target inflow
$T_i$ for each open edge:

1. **At every barotropic substep, after Flather's fill,** the transport on edge $i$ is shifted uniformly over its wet
   faces so that its integral equals the target (Oceananigans' `target_transport`):

   $$V_b \leftarrow V_b - \frac{\sum V_b\,\Delta x - s_i T_i}{L_i}, \qquad L_i = \sum_{\text{wet}} \Delta x .$$

   The stored averaged transport is pinned the same way at the end of the step. Flather still sets the *shape* of the
   flow along the edge, and still lets outgoing waves leave, but the water through each edge is exactly $T_i$.
2. **GLORYS's depth mean on each edge** is shifted so that $\sum s_i\,\bar v^G H\,\Delta x = T_i$. The velocities that
   carry T and S (2.3) then integrate to the same inflow as the water.

With $T_{Pacific} = +0.8$ Sv and $T_{Atlantic} = -0.8$ Sv (the Bering Strait throughflow, in through the Pacific
edge and out through the Atlantic), the edges add no net volume; the domain's volume changes only through the
surface freshwater flux (sea-ice growth and melt, precipitation, rivers).

### 4.4 The fix in use: pin only the total inflow (`net_inflow`, `NetInflow`)

The per-edge pin held the volume but piled water up at the Pacific edge (section 5, run 2): the model at 1° cannot
carry 0.8 Sv north from its Pacific edge. NEMO's `ln_vol` and ROMS's `VolCons` constrain only the *total* inflow over
all open boundaries, and we now do the same. Given one target $T$ for both edges together:

1. **At every barotropic substep, after Flather's fill,** one uniform correction is applied to the wet faces of both
   edges so that their total inflow equals $T$:

   $$\delta = \frac{\sum_i s_i \sum V_b\,\Delta x - T}{\sum_i L_i}, \qquad V_b \leftarrow V_b - s_i\,\delta \quad \text{on edge } i .$$

   The stored averaged transport is pinned the same way at the end of the step. Flather still decides how the inflow
   is *split* between the Atlantic and Pacific edges, through their sea-level mismatches, and still radiates outgoing
   waves; only the total is fixed.
2. **GLORYS's depth mean** is shifted by one common velocity on both edges so that its total exterior inflow is $T$
   (as `zero_net_inflow`, with target $T$ instead of 0):

   $$\bar v^G \leftarrow \bar v^G - s_i\,\frac{Q^{ext} - T}{A_b} .$$

In Oceananigans this is `GravityWaveRadiationBoundaryCondition((U, η); target_transport = NetInflow(T))` on every open
side that shares the target; in NumericalEarth, `nested_ocean_model(...; net_inflow = T)`; in the pan-Arctic config,
`[boundaries] net_inflow_Sv`. With $T = 0$ the edges add no net volume and the domain's volume changes only through the
surface freshwater flux.

### 4.5 The target: follow GLORYS's volume (`follow_glorys_volume`)

With $T = 0$ the domain loses the surface freshwater (−6 cm over January, mostly sea-ice growth) while GLORYS's mean η
falls only 3.5 cm, with day-to-day swings of ±2 cm. The target is now set before every time step from GLORYS's
domain-mean η, $\bar\eta^G(t)$ (its daily means sampled onto the model grid and interpolated linearly in time):

$$T = -F + A\,\frac{(\bar\eta^G - \bar\eta^G_0) - (\bar\eta - \bar\eta_0)}{\tau} \;\Big[\, + A\,\frac{d\bar\eta^G}{dt} \,\Big],$$

with $A$ the wet area, $F$ the surface freshwater volume flux into the ocean, $\bar\eta$ the model's domain-mean η and
the subscript 0 the start of the run. The first term hands the surface freshwater back through the edges; the second
relaxes the model's mean η toward GLORYS's over $\tau$ (1 day by default). The bracketed tendency term is optional
(`follow_glorys_tendency`): with it the model follows GLORYS exactly, but GLORYS's slope is constant within a day and
jumps at each new daily mean, so the inflow jumps by up to 5 Sv once a day. Relaxation alone is smooth and lags GLORYS
by about $\tau$; it is the default. The pinned value is held in a `Ref` (`NetInflow(Ref(T))`) that the pin reads at
every barotropic step. The shift of GLORYS's depth mean (4.4, point 2) uses the target at each new GLORYS sample.

### 4.6 Comparison

| | Volume through the edges | Velocities carrying T, S | Where the correction acts |
|---|---|---|---|
| #676 (Flather) | set by the sea-level mismatch (3.1) | agree with the water (3.3) | — |
| + `zero_net_inflow` | as above, minus the net of GLORYS's term | agree with the water | GLORYS data, once per sample |
| + outside control | held, at the cost of ±100 Sv shifts | inconsistent by ±40 Sv | GLORYS data, once per step |
| + `edge_transports` | exactly $T_i$ per edge; edge η drifts (+1.2 m) | integrate to $T_i$ | inside every substep, per edge |
| + `net_inflow` | total exactly $T$; Flather splits it | agree with the water | inside every substep, total (NEMO, ROMS) |
| + GLORYS target | total follows GLORYS's volume (±1 cm) | agree with the water | as above; $T$ set every step |

## 5. Evidence and tests

- **Volume budget** (monitor, every time step): the change in domain volume is split into the part through the edges,
  the part across the surface, and a residual. With #676 and with `zero_net_inflow`, the residual stays within
  0.4 cm over 30 days: the budget closes. The ±10 cm swings come through the edges; the surface removes a steady 6 cm
  over January (sea-ice growth).
- **Unit tests** (not yet run in the test suites):
  - NumericalEarth `test/test_nested_ocean_model.jl`: with 3 Sv in through one side and 3 Sv out through another
    (`edge_transports`), each side carries its target to within 10⁻⁶ after three steps and the volume is unchanged;
    with `net_inflow` = 2 Sv, the two sides together carry 2 Sv and the volume rises by 2 Sv × the elapsed time.
  - Oceananigans `test/hydrostatic_free_surface/open_boundary_conditions.jl`: two sides sharing `NetInflow(0)` carry
    no net inflow while each still has flow; `NetInflow(Q/2)` raises mean η by $Q\,\Delta t / 2A$.
- **30-day coupled runs** (1°, January 1993, JRA55 winds, sea ice, pressure off):
  1. **plain #676** (stopped at day 26.5, out of GPU memory after a previous run in the same session): the budget
     closes (residual within 0.4 cm); net inflow through the edges swings between −14 and +11 Sv; mean η rises 12 cm
     above its start by day 16 and falls 14 cm below by day 26, against GLORYS's +2 and −4 cm. The edge velocities and
     the volume agree (3.3).
  2. **`edge_transports`, +0.8 Sv in through the Pacific edge and 0.8 Sv out through the Atlantic edge:** the net
     inflow is zero to 0.005 Sv all month; mean η follows the surface freshwater curve exactly (−6.0 cm at day 30
     against GLORYS's −3.5 cm), residual zero. **But the model's sea level at the Pacific edge drifts far above
     GLORYS's: +17 cm on day 1, +96 cm on day 20, +122 cm on day 30** (Atlantic edge −3 to −14 cm), and the
     largest speed reaches 2.1 m/s. With the inflow pinned, Flather no longer ties the edge sea level to GLORYS's,
     and 0.8 Sv is more than the 1° model carries northward from its Pacific edge (about 0.2 Sv in the plain run):
     water piles up between the edge and Bering Strait.
  3. **`net_inflow` = 0 (total only):** the net inflow is zero from the first step on; mean η follows the surface
     freshwater curve exactly (−6.0 cm at day 30, against GLORYS's −3.5 cm, and up to 4 cm from GLORYS in between);
     residual zero; the edge velocities and the volume agree. Flather now splits the flow: the Atlantic and Pacific
     inflows vary by ±5.7 Sv, with a mean of 0.24 Sv in through the Pacific edge. The edge sea levels stay within
     −20 to +24 cm of GLORYS's without drifting (−5 and +2 cm at day 30).
  4. **GLORYS target with its tendency** (τ = 2 days): mean η follows GLORYS's to within 0.01 cm all month; the edges
     carry GLORYS's volume change minus the surface freshwater (+2.5 cm over the month). The net inflow jumps by up
     to 5.6 Sv at each new daily mean (between −2.8 and +4.8 Sv).
  5. **GLORYS target by relaxation alone** (τ = 1 day, the default): mean η follows GLORYS's about a day behind,
     within ±1.1 cm (rms 0.45 cm, −0.2 cm at day 30); the net inflow is smooth (largest hourly change 0.27 Sv,
     between −2.4 and +3.2 Sv). Edge sea levels: Atlantic −13 to +16 cm, Pacific −10 to +29 cm, no drift.

  The size of the net inflow in runs 4 and 5 is what GLORYS's volume requires: over this domain, 1 cm of mean η per
  day is 2.4 Sv. It is small against the plain run's −14 to +11 Sv and against the gross flow through each edge.

  Common to runs 3–5, so not set by the volume control: the largest speed reaches 2.3 m/s (identical in all three);
  the shelf (H < 200 m) falls about 10 cm while the deep ocean (H > 1000 m) rises about 1 cm.

## 6. Still open

- **Move the GLORYS target into NumericalEarth.** It is a script in the pan-Arctic sandbox (`glorys_volume_target.jl`);
  the nested model could take it as an option, since it holds the parent.
- **Weighting of the pin.** The pin shifts the transport uniformly per unit length of edge, which is a larger velocity
  shift on shallow faces. ROMS's `VolCons` shifts the velocity uniformly (weighted by $H\,\Delta x$), as our shift of
  GLORYS's depth mean already does; making the pin match is a small change.
- **The largest speeds (2.3 m/s) and the shelf drain** (section 5): to be located and compared with GLORYS on the
  model grid. The monitor now records the grid point of the largest speed every hour (`plot_max_speed.jl`).
- **Atmospheric pressure.** GLORYS's η has no inverse-barometer response, so pressure is off. Turning it on requires
  adding the inverse barometer to GLORYS's η at the edges and in the initial state (as NEMO's `ln_apr_obc`).
- **Initial state.** GLORYS interpolated, not balanced: the first ~8 hours swing by ±20 Sv in every run. The
  interpolated velocities at the edge faces carry about 22 Sv of net inflow before the first step; the boundary
  conditions remove it in the first step, a sudden adjustment. A balanced initial state would start consistent with
  the target.
- **Edge vertical structure.** Compare the model's and GLORYS's velocity and temperature sections along the edges,
  in particular the Atlantic Water layer.

## Where the code is

- Flather and Chapman: Oceananigans `src/BoundaryConditions/open_boundary_schemes/gravity_wave_schemes.jl`.
- `NetInflow` and the total pin: Oceananigans
  `src/Models/HydrostaticFreeSurfaceModels/SplitExplicitFreeSurfaces/barotropic_targeted_transport.jl` and
  `src/BoundaryConditions/open_boundary_schemes/open_boundary_utils.jl`.
- Substep loop, the edge pin and the end-of-step refill: Oceananigans
  `src/Models/HydrostaticFreeSurfaceModels/SplitExplicitFreeSurfaces/step_split_explicit_free_surface.jl`
  (substeps lines 144–163; refill and re-pin lines 288–296).
- Edge conditions built from GLORYS: NumericalEarth `src/NestedModels/ocean_open_boundaries.jl`.
- GLORYS sampling, `zero_net_inflow`, `net_inflow` and the per-edge shift: NumericalEarth
  `src/NestedModels/ocean_boundary_exchanger.jl`.
- GLORYS target: PanArcticShenanigans `nested/glorys_volume_target.jl`.
- Monitor and volume budget: PanArcticShenanigans `nested/monitor.jl`; plots `nested/plot_budget.jl`,
  `nested/plot_max_speed.jl`.
