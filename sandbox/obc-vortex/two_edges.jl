# Two open edges: a toy Atlantic (the whole south side, 2000 km) and a toy Bering Strait (a 300 km wide, 100 km long
# strait in the north wall), with a 1 Sv throughflow that enters at the strait and leaves through the south side.
# Flat bottom (H = 500 m), the vortex of obc_vortex.jl on the south side, the throughflow ramped up over 2 days.
#
# The external (parent) data are geostrophic: along each edge V = (gH/f) ∂η/∂x, so the throughflow comes with an η
# drop of fQ/(gH) = 2.9 cm across each edge, high on the west wall and low on the east wall.
#
# Cases (the error is in the parent data; the truth is the base run):
#   :base             consistent parent data
#   :datum_both       the parent η 3 cm too high at both edges (one datum for the whole parent model)
#   :datum_south      3 cm too high at the south edge only (datums that differ between edges)
#   :parent_mismatch  the parent's net transports do not close: 1.05 Sv out at the south, 1 Sv in at the strait
#
# Methods (all Flather underneath, V = Vᵉˣᵗ ∓ c(ηᵇ − ηᵉˣᵗ), c = √(gH)):
#   :flather       nothing else
#   :edge_net      each edge's net transport set to the parent's net for that edge (M5b, per edge)
#   :total_zero    NEMO/ROMS: one uniform normal velocity on all open faces sets the total net inflow to the
#                  surface freshwater flux (zero here)
#   :hybrid        :edge_net, then the :total_zero correction (each edge follows the parent, the total closes our
#                  own volume budget)
#   :hybrid_H      as :hybrid, but the per-edge correction is a uniform velocity (transport change ∝ H, as NEMO's)
#                  rather than a uniform η shift (∝ √(gH))
#
# Variant `shallow`: a 100 m deep strait and the 50 m shelf of E4 along the west of the south edge, so the two
# per-edge weightings differ; methods :flather, :total_zero, :hybrid, :hybrid_H.
#
#   julia -t 8 --project=$OCEAN_HPC_ROOT two_edges.jl [shallow]          # results/two_edges_*.jls, table on screen
#   julia --project=$OCEAN_HPC_ROOT/viz two_edges_plot.jl [shallow]      # figures/two_edges[_shallow].png

include("obc_vortex.jl")

const P2 = merge(default_parameters(), (strait_x = (1400kilometers, 1700kilometers), strait_length = 100kilometers,
                                        Q = 1e6, Q_ramp = 2days, datum = 0.03, datum_ramp = 1day, stop = 20days,
                                        strait_depth = 500.0, shelf = false))
const P2shallow = merge(P2, (strait_depth = 100.0, shelf = true))

const variant = isempty(ARGS) ? "flat" : ARGS[1]
const two_edge_methods = variant == "shallow" ? (:flather, :total_zero, :hybrid, :hybrid_H) :
                                                (:flather, :edge_net, :total_zero, :hybrid)
const two_edge_cases   = (:base, :datum_both, :datum_south, :parent_mismatch)

#####
##### Grid: the basin plus a strait through the north wall
#####

in_strait(x, P) = P.strait_x[1] < x < P.strait_x[2]

function two_edge_grid(P)
    Nx = round(Int, P.L / P.Δ)
    Ly = P.L + P.strait_length
    Ny = round(Int, Ly / P.Δ)
    grid = RectilinearGrid(CPU(); size = (Nx, Ny, 1), x = (0, P.L), y = (0, Ly), z = (-P.H, 0),
                           halo = (4, 4, 2), topology = (Bounded, TwoSided{Open, Open}, Bounded))
    bottom(x, y) = y > P.L ? (in_strait(x, P) ? -P.strait_depth : 0.0) : -depth(x, P; shelf = P.shelf)
    P.shelf || return ImmersedBoundaryGrid(grid, GridFittedBottom(bottom))
    return ImmersedBoundaryGrid(grid, PartialCellBottom(bottom; minimum_fractional_cell_height = 0.05))
end

#####
##### Parent (external) data on both edges
#####

function two_edge_parameters(grid, P, method, case)
    x  = xnodes(grid, Center())
    W  = P.strait_x[2] - P.strait_x[1]
    Hb = [depth(xi, P; shelf = P.shelf) for xi in x]
    ΔηS = P.f * P.Q * P.L / (P.g * sum(Hb) * P.Δ)          # geostrophic η drop across each edge
    ΔηN = P.f * P.Q / (P.g * P.strait_depth)
    strait = [in_strait(xi, P) for xi in x]
    return (; Nx = length(x), Δ = P.Δ, L = P.L, W, g = P.g, strait,
            VS₀ = Hb .* [vortex_v(xi, 0.0, P) for xi in x],
            qS = Hb ./ (sum(Hb) * P.Δ),                       # the throughflow spread as a uniform velocity
            ηS₀ = [vortex_η(xi, 0.0, P) for xi in x],
            ηSᵗᶠ = ΔηS .* (0.5 .- x ./ P.L),
            ηNᵗᶠ = [s ? ΔηN * (0.5 - (xi - P.strait_x[1]) / W) : 0.0 for (s, xi) in zip(strait, x)],
            QS = case === :parent_mismatch ? 1.05P.Q : P.Q, QN = P.Q,
            δS = case in (:datum_both, :datum_south) ? P.datum : 0.0,
            δN = case === :datum_both ? P.datum : 0.0,
            τQ = P.Q_ramp, τδ = P.datum_ramp,
            edge = method in (:edge_net, :hybrid), edge_H = method === :hybrid_H,
            total = method in (:total_zero, :hybrid, :hybrid_H))
end

@inline ramp(t, τ) = min(t / τ, one(t))
# V is positive northward: out of the basin through the strait, into it through the south side
@inline Vext_south(i, t, p) = @inbounds p.VS₀[i] - ramp(t, p.τQ) * p.QS * p.qS[i]
@inline ηext_south(i, t, p) = @inbounds p.ηS₀[i] + ramp(t, p.τQ) * p.ηSᵗᶠ[i] + p.δS * ramp(t, p.τδ)
@inline Vext_north(i, t, p) = @inbounds ifelse(p.strait[i], -ramp(t, p.τQ) * p.QN / p.W, zero(t))
@inline ηext_north(i, t, p) = @inbounds ifelse(p.strait[i], ramp(t, p.τQ) * p.ηNᵗᶠ[i] + p.δN * ramp(t, p.τδ), zero(t))

# Sums along both edges: c-weighted mismatch and its weight (for the per-edge constant C), wet area (for the uniform
# velocity), and the Flather transports (south: into the basin; north: out of it).
@inline function edge_sums(grid, η, t, p)
    kη = grid.Nz + 1
    jN = grid.Ny + 1
    nS = dS = nN = dN = AS = AN = FS = FN = zero(grid)
    for i in 1:p.Nx
        HS = max(static_column_depthᶜᶠᵃ(i, 1, grid), zero(grid))
        HN = max(static_column_depthᶜᶠᵃ(i, jN, grid), zero(grid))
        cS = sqrt(p.g * HS)
        cN = sqrt(p.g * HN)
        mS = @inbounds (η[i, 0, kη] + η[i, 1, kη]) / 2 - ηext_south(i, t, p)
        mN = @inbounds (η[i, jN - 1, kη] + η[i, jN, kη]) / 2 - ηext_north(i, t, p)
        nS += cS * mS * p.Δ; dS += cS * p.Δ; AS += HS * p.Δ
        nN += cN * mN * p.Δ; dN += cN * p.Δ; AN += HN * p.Δ
        FS += ifelse(HS > 0, (Vext_south(i, t, p) - cS * mS) * p.Δ, zero(grid))
        FN += ifelse(HN > 0, (Vext_north(i, t, p) + cN * mN) * p.Δ, zero(grid))
    end
    return (; nS, dS, nN, dN, AS, AN, FS, FN)
end

# Per-edge corrections, either η shifts (CS, CN) or uniform velocities (uS, uN), then the uniform velocity w removed
# from the total inflow
@inline function corrections(grid, η, t, p)
    s  = edge_sums(grid, η, t, p)
    CS = ifelse(p.edge, s.nS / s.dS, zero(grid))
    CN = ifelse(p.edge, s.nN / s.dN, zero(grid))
    uS = ifelse(p.edge_H, -s.nS / s.AS, zero(grid))
    uN = ifelse(p.edge_H, s.nN / s.AN, zero(grid))
    inflow = (s.FS + CS * s.dS - uS * s.AS) - (s.FN - CN * s.dN - uN * s.AN)
    w  = ifelse(p.total, inflow / (s.AS + s.AN), zero(grid))
    return CS, CN, uS, uN, w
end

@inline function south_value(i, k, grid, clock, fields, p)
    CS, CN, uS, uN, w = corrections(grid, fields.η, clock.time, p)
    H = max(static_column_depthᶜᶠᵃ(i, 1, grid), zero(grid))
    return (Vext_south(i, clock.time, p) - H * (uS + w), ηext_south(i, clock.time, p) + CS)
end

@inline function north_value(i, k, grid, clock, fields, p)
    CS, CN, uS, uN, w = corrections(grid, fields.η, clock.time, p)
    H = max(static_column_depthᶜᶠᵃ(i, grid.Ny + 1, grid), zero(grid))
    return (Vext_north(i, clock.time, p) - H * (uN - w), ηext_north(i, clock.time, p) + CN)
end

#####
##### Run
#####

function run_two_edges(method, case; P = P2)
    grid = two_edge_grid(P)
    p = two_edge_parameters(grid, P, method, case)
    V = FieldBoundaryConditions(south = GravityWaveRadiationBoundaryCondition(south_value; discrete_form = true, parameters = p),
                                north = GravityWaveRadiationBoundaryCondition(north_value; discrete_form = true, parameters = p))
    vbc = NormalFlowBoundaryCondition(0; scheme = NormalRadiation(inflow_timescale = Inf, outflow_timescale = Inf))
    v = FieldBoundaryConditions(south = vbc, north = vbc)
    free_surface = SplitExplicitFreeSurface(grid; substeps = P.substeps, gravitational_acceleration = P.g)
    model = HydrostaticFreeSurfaceModel(grid; free_surface, coriolis = FPlane(f = P.f), momentum_advection = nothing,
                                        buoyancy = nothing, tracers = (), boundary_conditions = (; V, v))
    set!(model, u = (x, y, z) -> vortex_u(x, y, P), v = (x, y, z) -> vortex_v(x, y, P), η = (x, y, z) -> vortex_η(x, y, P))

    Nx, Ny = size(grid)[1:2]
    xc, yc = xnodes(grid, Center()), ynodes(grid, Center())
    wet = [yc[j] < P.L || in_strait(xc[i], P) for i in 1:Nx, j in 1:Ny]
    Vb = model.free_surface.barotropic_velocities.V
    η  = model.free_surface.displacement
    times = Float64[]; south = Float64[]; north = Float64[]; mean_η = Float64[]
    function record!(sim)
        Va = Array(interior(Vb))[:, :, 1]
        ηa = Array(interior(η))[:, :, 1]
        push!(times, sim.model.clock.time)
        push!(south, sum(Va[:, 1]) * P.Δ)              # into the basin
        push!(north, sum(Va[:, Ny + 1]) * P.Δ)         # out of the basin
        push!(mean_η, sum(ηa[wet]) / count(wet))
        return nothing
    end
    simulation = Simulation(model; Δt = P.Δt, stop_time = P.stop, verbose = false)
    add_callback!(simulation, record!, TimeInterval(1hour))
    run!(simulation)
    result = (; method, case, times, south, north, mean_η, η_end = Array(interior(η))[:, :, 1], wet)
    mkpath("results")
    serialize("results/two_edges_$(variant == "flat" ? "" : variant * "_")$(method)_$(case).jls", result)
    return result
end

const Pvariant = variant == "shallow" ? P2shallow : P2
R = Dict((m, c) => run_two_edges(m, c; P = Pvariant) for c in two_edge_cases for m in two_edge_methods)

last_day(r, a) = (n = findall(t -> t ≥ r.times[end] - 1day, r.times); sum(a[n]) / length(n))
println("\nTwo open edges ($variant), day 20 (last-day means). Parent: 1 Sv in at the strait, 1 Sv out at the south",
        " (1.05 Sv for parent_mismatch). The truth for every case is the base run.")
@printf("  %-16s %-11s %10s %10s %12s %14s\n", "case", "method", "south [Sv]", "strait [Sv]", "Δ⟨η⟩ [cm]",
        "rms vs base [cm]")
for c in two_edge_cases, m in two_edge_methods
    r, b = R[(m, c)], R[(m, :base)]
    Δ = r.η_end[r.wet] .- b.η_end[b.wet]
    @printf("  %-16s %-11s %+10.3f %+10.3f %+12.2f %14.2f\n", c, m, -last_day(r, r.south) / 1e6,
            last_day(r, r.north) / 1e6 * -1, 100 * (r.mean_η[end] - r.mean_η[1]), 100 * sqrt(sum(abs2, Δ) / length(Δ)))
end
println("  (south: net outflow, positive out of the basin; strait: net inflow, positive into the basin)")
println("\nPlot with:  julia --project=\$OCEAN_HPC_ROOT/viz two_edges_plot.jl $(variant == "flat" ? "" : variant)")
