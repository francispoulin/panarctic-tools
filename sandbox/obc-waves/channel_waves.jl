# Rungs 0–2 of the open-boundary validation in Oceananigans: an x-periodic channel (no variation along the open edge,
# so this is the edge-mean mode), linear rotating shallow water, open to the south, wall to the north.
#
#   rung 0: a steady balanced along-edge jet, with the external η correct and 2 cm too high (the datum)
#   rung 1: one southward Poincaré wave packet at ω/f = 1.3, 1.6, 2, 3, 5; measured reflection against
#           r = (cₚ − α)/(cₚ + α), cₚ = c ω/√(ω² − f²), α = c for Flather, 0 for Dirichlet
#
# The same tests in Python (sw1d.py and friends) give obc_1d_validation.png; this checks that Oceananigans'
# GravityWaveRadiation reproduces them.
#
#   julia -t 8 --project=$OCEAN_HPC_ROOT channel_waves.jl            # results/channel.jls, a table on screen
#   julia --project=$OCEAN_HPC_ROOT/viz channel_plot.jl              # figures/channel_validation.png

using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: GravityWaveRadiationBoundaryCondition, NormalRadiation
using Oceananigans.Grids: TwoSided, Open, Closed
using Printf, Serialization

const g  = 9.80665
const H  = 500.0
const f  = 1.4e-4
const c  = sqrt(g * H)
const Δ  = 10kilometers
const Nx = 4                 # periodic in x: nothing varies along the edge
const Δt = 10minutes

phase_speed(ω) = c * ω / sqrt(ω^2 - f^2)
theory_r(ω, α) = (phase_speed(ω) - α) / (phase_speed(ω) + α)

#####
##### Boundary conditions on the barotropic transport V at the south face
#####

@inline flather_value(i, k, grid, clock, fields, p) = (zero(grid), p.ηᵉˣᵗ)

struct ZeroTransport end                                   # Dirichlet with Vᵉˣᵗ = 0 (both call forms)
@inline (::ZeroTransport)(i, k, grid) = zero(grid)
@inline (::ZeroTransport)(i, k, grid, clock, fields) = zero(grid)

function boundary_conditions(method, ηᵉˣᵗ)
    V = method === :flather ?
        GravityWaveRadiationBoundaryCondition(flather_value; discrete_form = true, parameters = (; ηᵉˣᵗ)) :
        NormalFlowBoundaryCondition(ZeroTransport())
    v = NormalFlowBoundaryCondition(0; scheme = NormalRadiation(inflow_timescale = Inf, outflow_timescale = Inf))
    return (V = FieldBoundaryConditions(south = V), v = FieldBoundaryConditions(south = v))
end

function build_model(method; Ly, Lext = 0.0, reference = false, ηᵉˣᵗ = 0.0)
    y₁ = reference ? -Lext : 0.0
    Ny = round(Int, (Ly - y₁) / Δ)
    TY = reference ? Bounded : TwoSided{Open, Closed}
    grid = RectilinearGrid(CPU(); size = (Nx, Ny, 1), x = (0, Nx * Δ), y = (y₁, Ly), z = (-H, 0),
                           halo = (4, 4, 1), topology = (Periodic, TY, Bounded))
    bcs = reference ? NamedTuple() : boundary_conditions(method, ηᵉˣᵗ)
    free_surface = SplitExplicitFreeSurface(grid; substeps = 12, gravitational_acceleration = g)
    return HydrostaticFreeSurfaceModel(grid; free_surface, coriolis = FPlane(; f), momentum_advection = nothing,
                                       buoyancy = nothing, tracers = (), boundary_conditions = bcs)
end

η_profile(model) = Array(interior(model.free_surface.displacement))[1, :, 1]

function run_record!(model, T, ys)
    yc = ynodes(model.grid, Center())
    js = argmin(abs.(yc .- ys))
    η = model.free_surface.displacement
    rec = Float64[]
    record!(sim) = push!(rec, interior(η)[1, js, 1])
    simulation = Simulation(model; Δt, stop_time = T, verbose = false)
    add_callback!(simulation, record!, IterationInterval(1))
    run!(simulation)
    return rec
end

#####
##### Rung 1: one southward Poincaré packet
#####

function packet(ω; A = 0.01)
    l = sqrt(ω^2 - f^2) / c
    σ = 1.5 * 2π / l
    return (; ω, l, Y = ω / l, A, σ, yp = 4σ + 500kilometers)
end
envelope(p, y) = p.A * exp(-(y - p.yp)^2 / (2p.σ^2))
packet_η(p, y) = envelope(p, y) * cos(p.l * (y - p.yp))
packet_u(p, y) = -(f * p.Y / (p.ω * H)) * envelope(p, y) * sin(p.l * (y - p.yp))
packet_v(p, y) = -(p.Y / H) * envelope(p, y) * cos(p.l * (y - p.yp))          # southward: V = −cₚ η

function reflection(ω; ys = 500kilometers)
    p  = packet(ω)
    cg = c^2 * p.l / ω
    T  = ceil(1.15 * (p.yp + 4p.σ + 2ys) / cg / Δt) * Δt
    Lext = ceil(c * T / Δ) * Δ
    Ly   = ceil((p.yp + 4p.σ + c * T / 2) / Δ) * Δ
    out = Dict{Symbol, Any}()
    for method in (:reference, :flather, :dirichlet)
        model = build_model(method; Ly, Lext, reference = method === :reference)
        set!(model, u = (x, y, z) -> packet_u(p, y), v = (x, y, z) -> packet_v(p, y), η = (x, y, z) -> packet_η(p, y))
        out[method] = run_record!(model, T, ys)
    end
    r(m) = sqrt(sum(abs2, out[m] .- out[:reference]) / sum(abs2, out[:reference]))
    return (; ω, flather = r(:flather), dirichlet = r(:dirichlet), series = out)
end

#####
##### Rung 0: steady balanced jet, external η correct or 2 cm high
#####

const jet_A, jet_y, jet_w = 0.1, 600kilometers, 400kilometers
jet_η(y) = jet_A * exp(-(y - jet_y)^2 / (2jet_w^2))
jet_u(y) = (g / f) * jet_η(y) * (y - jet_y) / jet_w^2                         # u = −(g/f) ∂η/∂y

function steady(method, δ; Ly = 3000kilometers, T = 5days)
    model = build_model(method; Ly, ηᵉˣᵗ = jet_η(0.0) + δ)
    set!(model, u = (x, y, z) -> jet_u(y), η = (x, y, z) -> jet_η(y))
    run!(Simulation(model; Δt, stop_time = T, verbose = false))
    y  = ynodes(model.grid, Center())
    η  = η_profile(model)
    u  = Array(interior(model.velocities.u))[1, :, 1]
    return (; y, η, truth = jet_η.(y), err = sqrt(sum(abs2, η .- jet_η.(y)) / length(y)),
            jet_Sv = H * sum(u .- jet_u.(y)) * Δ / 1e6)
end

#####
##### Run everything
#####

results = Dict{Symbol, Any}()
println("Rung 0: steady balanced jet, day 5")
for δ in (0.0, 0.02), method in (:flather, :dirichlet)
    s = steady(method, δ)
    results[Symbol(method, "_", round(Int, 100δ), "cm")] = s
    @printf("  datum %.0f cm  %-9s  rms η error %.3f cm   spurious along-edge transport %+.3f Sv\n",
            100δ, method, 100s.err, s.jet_Sv)
end
@printf("  expected for Flather with a 2 cm datum: gHδ/f = %.3f Sv\n", g * H * 0.02 / f / 1e6)

println("\nRung 1: reflection of one Poincaré packet (η), measured and theory")
@printf("  %6s  %16s  %16s\n", "ω/f", "Flather (theory)", "Dirichlet")
packets = []
for r in (1.3, 1.6, 2.0, 3.0, 5.0)
    res = reflection(r * f)
    push!(packets, res)
    @printf("  %6.1f  %7.3f (%5.3f)  %16.3f\n", r, res.flather, abs(theory_r(r * f, c)), res.dirichlet)
end
results[:packets] = [(; ω = p.ω, flather = p.flather, dirichlet = p.dirichlet) for p in packets]

mkpath("results")
serialize("results/channel.jls", results)
println("\nsaved results/channel.jls; plot with  julia --project=\$OCEAN_HPC_ROOT/viz channel_plot.jl")
