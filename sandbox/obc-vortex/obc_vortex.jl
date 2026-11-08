# Open boundary conditions for the barotropic transport: a Gaussian vortex in a basin with one open side.
#
# Linear rotating shallow water as a HydrostaticFreeSurfaceModel with Nz = 1, no advection, no buoyancy, an f-plane and a
# SplitExplicitFreeSurface (so the open-boundary machinery is the one used in the pan-Arctic model). Walls on the west,
# north and east; the south side (y = 0) is open. A geostrophic Gaussian vortex centred on the open side crosses it
# inwards on one half and outwards on the other, with zero net transport. See balanced_obc_notes.pdf.
#
# Boundary conditions on the barotropic transport V at the south faces (method symbols):
#   :flather         M1  V = Vᵉˣᵗ − c (ηᵇ − ηᵉˣᵗ)
#   :dirichlet       M2  V = Vᵉˣᵗ
#   :davies          M3  M2 plus relaxation of u, v to the external state in a rim
#   sponge = true adds the same rim to any method (Flather + sponge is NEMO's default set-up)
#   :flather_target  M4  M1 with target_transport = net external transport (uniform shift per unit length)
#   :balanced_net    M5b V = Vᵉˣᵗ − c (ηᵇ − ηᵉˣᵗ − C), C from the net constraint (√(gH)-weighted)
#   :balanced_pv     M5a same, C from the PV inversion with the volume constraint (updated every model step)
#   :balanced_filtered M6 same, C a slow average (timescale T_C) of the √(gH)-weighted mean edge mismatch
#   :rotating_flather M7 V = Vᵉˣᵗ − c (ηᵇ − ηᵉˣᵗ) + S/2,  dS/dt = −f H (uᵇ − uᵉˣᵗ) at each edge point: Flather with
#                        the rotational correction to the wave admittance (sandbox/obc-waves); S also absorbs a datum
#
# Experiments: :steady (E0), :adjustment (E1), :offset (E2), :offset_ramp (E2b: the offset grows over δ_ramp),
# :offset_slow (E2c: the offset grows over δ_ramp_slow, run for 20 days),
# :net_inflow (E3); keyword shelf = true for E4.
#
# Usage:
#   include("obc_vortex.jl")
#   result = run_case(:flather, :adjustment)          # returns a NamedTuple, also serialised to results/
#   reference = run_reference(:adjustment)

using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: GravityWaveRadiationBoundaryCondition, NormalRadiation
using Oceananigans.Grids: static_column_depthᶜᶠᵃ, TwoSided, Open, Closed
using Oceananigans.Operators: Δxᶜᶠᶜ
using Printf, Serialization

include("pv_inversion.jl")

#####
##### Parameters
#####

const g_earth = 9.80665

default_parameters() = (
    L   = 2000kilometers,       # basin size
    Δ   = 20kilometers,         # grid spacing
    H   = 500.0,                # depth (m): c = 70 m/s, L_d = 500 km
    f   = 1.4e-4,               # Coriolis parameter (≈ 70°N)
    g   = g_earth,
    Δt  = 20minutes,
    substeps = 30,
    # vortex (anticyclone for η₀ > 0), centred on the open side
    η₀  = 0.1, R = 150kilometers, x₀ = 1000kilometers, y₀ = 0.0,
    # E1: unbalanced bump (η only, no velocity)
    bump_A = 0.2, bump_R = 150kilometers, bump_x = 1000kilometers, bump_y = 1500kilometers,
    # E2: datum offset of the external free surface
    δ   = 0.03,
    δ_ramp = 1day,              # :offset_ramp: the offset grows linearly over this time
    δ_ramp_slow = 10days,       # :offset_slow: the same, slowly (a realistic drift of the datum)
    # E3: slowly varying net inflow Q(t) = Q₀ sin(2πt/T_Q), spread uniformly along the open side
    Q₀  = 1e6, T_Q = 10days,
    # M3: Davies rim
    rim_width = 200kilometers, rim_rate = 1 / 2hours,
    # M6: averaging timescale for C
    T_C = 2days,
    # E4: shelf along the west part of the basin
    shelf_depth = 50.0, shelf_x = 300kilometers, shelf_w = 50kilometers,
    # reference domain: extended this far to the south (walls everywhere)
    south_extension = 8000kilometers,
)

#####
##### Vortex, bump, shelf
#####

vortex_η(x, y, P) = P.η₀ * exp(-((x - P.x₀)^2 + (y - P.y₀)^2) / P.R^2)
vortex_u(x, y, P) =  (P.g / P.f) * vortex_η(x, y, P) * 2 * (y - P.y₀) / P.R^2     # u = −(g/f) ∂η/∂y
vortex_v(x, y, P) = -(P.g / P.f) * vortex_η(x, y, P) * 2 * (x - P.x₀) / P.R^2     # v =  (g/f) ∂η/∂x
bump_η(x, y, P)   = P.bump_A * exp(-((x - P.bump_x)^2 + (y - P.bump_y)^2) / P.bump_R^2)

depth(x, P; shelf = false) = shelf ? P.H - (P.H - P.shelf_depth) * (1 - tanh((x - P.shelf_x) / P.shelf_w)) / 2 : P.H

# E3: net inflow and the volume it has brought in
net_inflow(t, P, experiment)    = experiment === :net_inflow ? P.Q₀ * sin(2π * t / P.T_Q) : 0.0
inflow_volume(t, P, experiment) = experiment === :net_inflow ? P.Q₀ * P.T_Q / 2π * (1 - cos(2π * t / P.T_Q)) : 0.0

#####
##### Grid
#####

# Oceananigans ≥ 0.117: open-boundary schemes need an `Open` side, so the basin is open to the south and closed to
# the north (`TwoSided{Open, Closed}` in y). The reference basin is closed on all sides.
function build_grid(P; reference = false, shelf = false)
    Nx = round(Int, P.L / P.Δ)
    y₁ = reference ? -P.south_extension : 0.0
    Ny = round(Int, (P.L - y₁) / P.Δ)
    TY = reference ? Bounded : TwoSided{Open, Closed}
    grid = RectilinearGrid(CPU(); size = (Nx, Ny, 1), x = (0, P.L), y = (y₁, P.L), z = (-P.H, 0),
                           halo = (4, 4, shelf ? 2 : 1), topology = (Bounded, TY, Bounded))
    shelf || return grid
    bottom(x, y) = -depth(x, P; shelf = true)
    return ImmersedBoundaryGrid(grid, PartialCellBottom(bottom; minimum_fractional_cell_height = 0.05))
end

#####
##### External data on the south faces
#####
#
# p.V₀[i], p.η₀[i]: steady external transport per unit width (m²/s) and free surface at the south faces.
# Time dependence (E3) and the datum offset (E2) are added in the functions below.

function external_parameters(grid, P, experiment; shelf = false)
    Nx = size(grid, 1)
    x  = xnodes(grid, Center())
    Hb = [depth(xi, P; shelf) for xi in x]
    V₀ = Hb .* [vortex_v(xi, 0.0, P) for xi in x]
    η₀ = [vortex_η(xi, 0.0, P) for xi in x]
    δ  = experiment in (:offset, :offset_ramp, :offset_slow) ? P.δ : 0.0
    τδ = experiment === :offset_ramp ? P.δ_ramp : experiment === :offset_slow ? P.δ_ramp_slow : 0.0
    E3 = experiment === :net_inflow
    uᵉˣᵗ = [vortex_u(xi, P.Δ / 2, P) for xi in x]          # along-edge velocity at the first row of cells (M7)
    return (; V₀, η₀, δ, τδ, E3, Q₀ = P.Q₀, ω_Q = 2π / P.T_Q, Lx = P.L, A = P.L^2, g = P.g, C = [0.0],
            uᵉˣᵗ, S = zeros(Nx), Hb)
end

@inline external_Q(t, p)   = ifelse(p.E3, p.Q₀ * sin(p.ω_Q * t), zero(t))
@inline external_vol(t, p) = ifelse(p.E3, p.Q₀ / p.ω_Q * (1 - cos(p.ω_Q * t)), zero(t))
@inline external_V(i, t, p) = @inbounds p.V₀[i] + external_Q(t, p) / p.Lx
@inline offset_ramp(t, p)  = ifelse(p.τδ > 0, min(t / p.τδ, one(t)), one(t))
@inline external_η(i, t, p) = @inbounds p.η₀[i] + p.δ * offset_ramp(t, p) + external_vol(t, p) / p.A

# M1, M4: (Vᵉˣᵗ, ηᵉˣᵗ)
@inline flather_value(i, k, grid, clock, fields, p) = (external_V(i, clock.time, p), external_η(i, clock.time, p))

# M7: Flather with Vᵉˣᵗ + S/2; S is stepped in a callback (see run_case)
@inline rotating_flather_value(i, k, grid, clock, fields, p) =
    (@inbounds(external_V(i, clock.time, p) + p.S[i] / 2), external_η(i, clock.time, p))

# M2, M3: Vᵉˣᵗ. A plain normal-flow condition is also filled without a clock (at initialisation), which a
# DiscreteBoundaryFunction does not accept, so use a small callable that handles both forms.
struct DirichletTransport{P}
    p :: P
end
@inline (d::DirichletTransport)(i, k, grid) = external_V(i, zero(grid), d.p)
@inline (d::DirichletTransport)(i, k, grid, clock, fields) = external_V(i, clock.time, d.p)

# M5b: Flather about ηᵉˣᵗ + C, with C chosen so that the net transport equals the external net transport.
@inline function balanced_net_value(i, k, grid, clock, fields, p)
    t  = clock.time
    η  = fields.η
    kη = grid.Nz + 1
    num = zero(grid)
    den = zero(grid)
    for i′ in 1:grid.Nx
        c  = sqrt(p.g * max(static_column_depthᶜᶠᵃ(i′, 1, grid), zero(grid)))
        Δx = Δxᶜᶠᶜ(i′, 1, 1, grid)
        ηᵇ = (η[i′, 0, kη] + η[i′, 1, kη]) / 2
        num += c * (ηᵇ - external_η(i′, t, p)) * Δx
        den += c * Δx
    end
    C = num / den
    return (external_V(i, t, p), external_η(i, t, p) + C)
end

# M5a, M6: Flather about ηᵉˣᵗ + C, with C = p.C[1] updated by a callback every model step
# (M5a: from the PV inversion; M6: slow average of the √(gH)-weighted mean edge mismatch).
@inline balanced_pv_value(i, k, grid, clock, fields, p) =
    (external_V(i, clock.time, p), external_η(i, clock.time, p) + @inbounds(p.C[1]))

#####
##### Model
#####

function boundary_conditions(method, grid, p)
    V_south = if method === :flather
        GravityWaveRadiationBoundaryCondition(flather_value; discrete_form = true, parameters = p)
    elseif method === :flather_target
        Q = sum(p.V₀) * P_Δx(grid)
        GravityWaveRadiationBoundaryCondition(flather_value; discrete_form = true, parameters = p, target_transport = Q)
    elseif method === :balanced_net
        GravityWaveRadiationBoundaryCondition(balanced_net_value; discrete_form = true, parameters = p)
    elseif method === :rotating_flather
        GravityWaveRadiationBoundaryCondition(rotating_flather_value; discrete_form = true, parameters = p)
    elseif method in (:balanced_pv, :balanced_filtered)
        GravityWaveRadiationBoundaryCondition(balanced_pv_value; discrete_form = true, parameters = p)
    elseif method in (:dirichlet, :davies)
        NormalFlowBoundaryCondition(DirichletTransport(p))
    else
        error("unknown method $method")
    end
    # The three-dimensional v on the open faces: with Nz = 1 the barotropic correction sets it to V / H.
    v_south = NormalFlowBoundaryCondition(0; scheme = NormalRadiation(inflow_timescale = Inf, outflow_timescale = Inf))
    return (V = FieldBoundaryConditions(south = V_south), v = FieldBoundaryConditions(south = v_south))
end

P_Δx(grid) = grid isa ImmersedBoundaryGrid ? grid.underlying_grid.Δxᶜᵃᵃ : grid.Δxᶜᵃᵃ

function davies_forcing(P, experiment)
    mask(x, y, z) = y < P.rim_width ? (1 - y / P.rim_width)^2 : 0.0
    u_target(x, y, z, t) = vortex_u(x, y, P)
    v_target(x, y, z, t) = vortex_v(x, y, P) + net_inflow(t, P, experiment) / (P.L * depth(x, P))
    return (u = Relaxation(rate = P.rim_rate, mask = mask, target = u_target),
            v = Relaxation(rate = P.rim_rate, mask = mask, target = v_target))
end

function build_model(method, experiment; P = default_parameters(), shelf = false, reference = false, sponge = false)
    grid = build_grid(P; reference, shelf)
    p = external_parameters(grid, P, experiment; shelf)
    boundary_conditions_ = reference ? NamedTuple() : boundary_conditions(method, grid, p)
    forcing = (!reference && (method === :davies || sponge)) ? davies_forcing(P, experiment) : NamedTuple()
    free_surface = SplitExplicitFreeSurface(grid; substeps = P.substeps, gravitational_acceleration = P.g)
    model = HydrostaticFreeSurfaceModel(grid; free_surface, coriolis = FPlane(f = P.f),
                                        momentum_advection = nothing, buoyancy = nothing, tracers = (),
                                        boundary_conditions = boundary_conditions_, forcing)
    return model, p
end

function set_initial_state!(model, experiment, P)
    ηᵢ(x, y, z) = vortex_η(x, y, P) + (experiment === :adjustment ? bump_η(x, y, P) : 0.0)
    set!(model, u = (x, y, z) -> vortex_u(x, y, P), v = (x, y, z) -> vortex_v(x, y, P), η = ηᵢ)
    return model
end

#####
##### Run
#####

# Arrays of the interior state on the CPU: u (Nx+1, Ny), v (Nx, Ny+1), η (Nx, Ny)
state(model) = (u = Array(interior(model.velocities.u))[:, :, 1],
                v = Array(interior(model.velocities.v))[:, :, 1],
                η = Array(interior(model.free_surface.displacement))[:, :, 1])

function inverter_for(model, P; shelf = false, reference = false)
    grid = model.grid
    Nx, Ny = size(grid)[1:2]
    x = xnodes(grid, Center())
    H = [depth(x[i], P; shelf) for i in 1:Nx, j in 1:Ny]
    return PVInverter(Nx, Ny, P.Δ, P.Δ, P.f, P.g, H; open_south = !reference)
end

# The M5b constant evaluated from the model state (for the record; the boundary condition computes its own)
function net_constraint_C(model, p)
    grid = model.grid
    η  = model.free_surface.displacement
    kη = grid.Nz + 1
    t  = model.clock.time
    num = 0.0; den = 0.0
    for i in 1:grid.Nx
        c  = sqrt(p.g * max(static_column_depthᶜᶠᵃ(i, 1, grid), 0.0))
        ηᵇ = η[i, 1, kη]       # the first interior row: the halo row is still empty before the first time step
        num += c * (ηᵇ - external_η(i, t, p)); den += c
    end
    return num / den
end

stop_time_for(experiment) = experiment === :net_inflow ? 30days : experiment === :offset_slow ? 20days : 2days

"""
    run_case(method, experiment; P, shelf = false, stop_time, output_interval = 1hour, dir = "results")

Run one method on one experiment. Records hourly snapshots of (u, v, η), the net transport through the open side,
the domain-mean free surface and, for the balanced methods, the constant C. Returns a NamedTuple and serialises it
to `dir/<method>_<experiment>[_shelf].jls`.
"""
function run_case(method, experiment; P = default_parameters(), shelf = false, reference = false,
                  stop_time = stop_time_for(experiment), output_interval = 1hour, dir = "results", verbose = true,
                  suffix = "", sponge = false)
    model, p = build_model(method, experiment; P, shelf, reference, sponge)
    set_initial_state!(model, experiment, P)
    simulation = Simulation(model; Δt = P.Δt, stop_time, verbose = false)

    inv = inverter_for(model, P; shelf, reference)
    V = model.free_surface.barotropic_velocities.V
    Δx = P.Δ
    v_ext_now(t) = [external_V(i, t, p) / inv.H[i, 1] for i in 1:inv.Nx]

    # M5a: update C from the PV inversion every model step (lagged by one step)
    if method === :balanced_pv && !reference
        function update_C!(sim)
            s = state(sim.model)
            t = sim.model.clock.time
            B = balanced_state(inv, s.u, s.v, s.η; v_ext = v_ext_now(t))
            ηB_south = (P.f / P.g) .* B.sS                       # balanced η on the south faces
            ηext = [external_η(i, t, p) for i in 1:inv.Nx]
            p.C[1] = sum(ηB_south .- ηext) / inv.Nx
            return nothing
        end
        update_C!(simulation)
        add_callback!(simulation, update_C!, IterationInterval(1))
    end

    # M6: C relaxes towards the √(gH)-weighted mean edge mismatch m(t) over T_C (implicit exponential average).
    # At t = 0 any mismatch is taken to be a datum offset: C(0) = m(0).
    if method === :balanced_filtered && !reference
        p.C[1] = net_constraint_C(model, p)
        function filter_C!(sim)
            α = sim.Δt / P.T_C
            p.C[1] = (p.C[1] + α * net_constraint_C(sim.model, p)) / (1 + α)
            return nothing
        end
        add_callback!(simulation, filter_C!, IterationInterval(1))
    end

    # M7: dS/dt = −f H (uᵇ − uᵉˣᵗ), uᵇ the along-edge velocity at the centres of the first row of cells
    if method === :rotating_flather && !reference
        function step_S!(sim)
            u = sim.model.velocities.u
            for i in 1:length(p.S)
                uᵇ = (u[i, 1, 1] + u[i+1, 1, 1]) / 2
                p.S[i] -= sim.Δt * P.f * p.Hb[i] * (uᵇ - p.uᵉˣᵗ[i])
            end
            return nothing
        end
        add_callback!(simulation, step_S!, IterationInterval(1))
    end

    times = Float64[]; net = Float64[]; net_ext = Float64[]; mean_η = Float64[]; Cs = Float64[]
    snapshots = NamedTuple[]
    function record!(sim)
        t = sim.model.clock.time
        s = state(sim.model)
        push!(times, t)
        push!(net, reference ? 0.0 : sum(Array(interior(V))[:, 1, 1]) * Δx)
        push!(net_ext, sum(external_V(i, t, p) for i in 1:inv.Nx) * Δx)
        push!(mean_η, sum(s.η) / length(s.η))
        push!(Cs, method === :balanced_net ? net_constraint_C(sim.model, p) : p.C[1])
        push!(snapshots, map(a -> Float32.(a), s))
        verbose && @printf("  %-14s %-10s t = %6.2f d  net = %+7.3f Sv (ext %+7.3f)  ⟨η⟩ = %+.4f m\n",
                           reference ? "reference" : method, experiment, t / days, net[end] / 1e6, net_ext[end] / 1e6, mean_η[end])
        return nothing
    end
    add_callback!(simulation, record!, TimeInterval(output_interval))   # also fires at t = 0

    wall = @elapsed run!(simulation)

    name = (reference ? "reference" : String(method)) * suffix
    result = (; method = name, experiment = String(experiment), shelf, P, times, net, net_ext, mean_η,
              C = Cs, snapshots, wall)
    mkpath(dir)
    path = joinpath(dir, "$(name)_$(experiment)$(shelf ? "_shelf" : "").jls")
    serialize(path, result)
    verbose && @printf("  saved %s (%.1f s)\n", path, wall)
    return result
end

run_reference(experiment; kw...) = run_case(:none, experiment; reference = true, kw...)

const all_methods = (:flather, :dirichlet, :davies, :flather_target, :balanced_net, :balanced_pv, :balanced_filtered,
                     :rotating_flather)
