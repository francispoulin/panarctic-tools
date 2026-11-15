# Stratified flow over a ridge, with an energy budget. Without forcing, the total energy
#     E = Σ ½ (u² + v²) dV  +  Σ (−b z) dV  +  Σ ½ g η² dA          (per unit ρ₀)
# can only fall (WENO dissipates; its numerical mixing raises the potential energy a little). A partial-cell bottom
# that feeds energy into the flow shows up as E or the bottom-cell kinetic energy rising, and as fast flow in the
# bottom cells.
#
# Channel periodic in x, 32 cells of 15 km; a plateau of 8 cells at 1863 m in a 2837 m basin, so the flow crosses two
# 1000 m cliffs (as at the (306, 170) point of the 2005 run). Exponential levels to 5000 m, z-star, the GLORYS T and S
# profile with a linear equation of state (so that E is the exact energy of the equations). Initially a barotropic
# current of 0.1 m/s in the basin (constant transport across the ridge, 0.15 m/s over it); with rotation (f-plane at
# 72°N) a free-surface tilt across the channel balances the basin flow geostrophically.
#
# Cases: partial vs whole cells × 20 vs 50 levels × f = 0 vs 72°N, ten days each.
#
#   julia --project=<PanArcticShenanigans>/nested flow_over_ridge.jl [profile.csv] [days]

using Oceananigans, Oceananigans.Units, Printf
using Oceananigans.Operators: volume, Azᶜᶜᶜ
using Oceananigans.Grids: peripheral_node, inactive_node, znode
using Oceananigans.BuoyancyFormulations: buoyancy_perturbationᶜᶜᶜ

profile_file = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(@__DIR__, "glorys_profile_0E_72N.csv")
run_days     = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 10.0

rows = [parse.(Float64, split(line, ',')) for line in readlines(profile_file)[2:end] if !isempty(line)]
profile_depth = first.(rows)
profile_T     = getindex.(rows, 2)
profile_S     = getindex.(rows, 3)

function interpolate_profile(values, z)
    d = clamp(-z, first(profile_depth), last(profile_depth))
    k = clamp(searchsortedlast(profile_depth, d), 1, length(profile_depth) - 1)
    α = (d - profile_depth[k]) / (profile_depth[k + 1] - profile_depth[k])
    return (1 - α) * values[k] + α * values[k + 1]
end

Tᵢ(x, y, z) = interpolate_profile(profile_T, z)
Sᵢ(x, y, z) = interpolate_profile(profile_S, z)

Nx, Ny  = 32, 8
spacing = 15kilometers
Lx, Ly  = Nx * spacing, Ny * spacing
depth   = 5000
basin_depth, plateau_depth = 2837, 1863
on_plateau(x) = 12spacing ≤ x < 20spacing
bottom(x, y) = on_plateau(x) ? -plateau_depth : -basin_depth

U = 0.1
g = Oceananigans.defaults.gravitational_acceleration
uᵢ(x, y, z) = U * basin_depth / (on_plateau(x) ? plateau_depth : basin_depth)

function ridge_model(bottom_kind, Nz, latitude)
    z = ExponentialDiscretization(Nz, -depth, 0; scale = depth / 4, mutable = true)
    underlying_grid = RectilinearGrid(CPU(); size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z,
                                      halo = (7, 7, 7), topology = (Periodic, Bounded, Bounded))
    Bottom = bottom_kind == :partial ? PartialCellBottom : GridFittedBottom
    grid = ImmersedBoundaryGrid(underlying_grid, Bottom(bottom))

    coriolis = isnothing(latitude) ? nothing : FPlane(; latitude)
    f = isnothing(latitude) ? 0.0 : coriolis.f

    model = HydrostaticFreeSurfaceModel(grid;
                                        momentum_advection = WENOVectorInvariant(),
                                        tracer_advection = WENO(order = 5),
                                        buoyancy = SeawaterBuoyancy(equation_of_state = LinearEquationOfState(thermal_expansion = 2e-4, haline_contraction = 7.8e-4)),
                                        coriolis,
                                        free_surface = SplitExplicitFreeSurface(grid; substeps = 30),
                                        tracers = (:T, :S))

    set!(model; η = (x, y, z...) -> -f * U * (y - Ly / 2) / g)
    set!(model; u = uᵢ, T = Tᵢ, S = Sᵢ)
    return model
end

bottom_face(i, j, k, grid) = inactive_node(i - 1, j, k - 1, grid, Center(), Center(), Center()) |
                             inactive_node(i,     j, k - 1, grid, Center(), Center(), Center())

function energies(model)
    grid = model.grid
    Nz = size(grid, 3)
    u, v, η = model.velocities.u, model.velocities.v, model.free_surface.displacement
    kinetic = bottom_kinetic = potential = surface = 0.0
    fastest = (0.0, 0, 0, false)
    for i in 1:Nx, j in 1:Ny, k in 1:Nz
        if !peripheral_node(i, j, k, grid, Face(), Center(), Center())
            e = u[i, j, k]^2 / 2 * volume(i, j, k, grid, Face(), Center(), Center())
            kinetic += e
            at_bottom = bottom_face(i, j, k, grid)
            at_bottom && (bottom_kinetic += e)
            abs(u[i, j, k]) > abs(fastest[1]) && (fastest = (u[i, j, k], i, k, at_bottom))
        end
        peripheral_node(i, j, k, grid, Center(), Face(), Center()) ||
            (kinetic += v[i, j, k]^2 / 2 * volume(i, j, k, grid, Center(), Face(), Center()))
        inactive_node(i, j, k, grid, Center(), Center(), Center()) ||
            (potential -= buoyancy_perturbationᶜᶜᶜ(i, j, k, grid, model.buoyancy.formulation, model.tracers) *
                          znode(i, j, k, grid, Center(), Center(), Center()) * volume(i, j, k, grid, Center(), Center(), Center()))
    end
    for i in 1:Nx, j in 1:Ny
        surface += g * η[i, j, Nz + 1]^2 / 2 * Azᶜᶜᶜ(i, j, Nz + 1, grid)
    end
    return (; kinetic, bottom_kinetic, potential, surface, fastest)
end

println("Energy changes relative to the initial kinetic energy KE₀; max|u| in m/s at (i, k), * = a bottom face.")
println("The cliffs are the faces i = 13 and i = 21.\n")

for latitude in (nothing, 72), Nz in (20, 50), bottom_kind in (:partial, :whole)
    model = ridge_model(bottom_kind, Nz, latitude)
    E₀ = energies(model)
    KE₀ = E₀.kinetic
    total₀ = E₀.kinetic + E₀.potential + E₀.surface
    label = @sprintf("%-7s cells, Nz = %2d, %s", bottom_kind, Nz, isnothing(latitude) ? "f = 0   " : "f(72°N)")
    println(label)
    @printf("  %5s %9s %9s %9s %9s %9s %16s\n", "day", "ΔE", "ΔKE", "ΔKE_bot", "ΔPE", "ΔE_η", "max|u| (i, k)")
    Δt = 10minutes
    for n in 1:ceil(Int, run_days)
        for _ in 1:round(Int, Oceananigans.Units.day / Δt)
            time_step!(model, Δt)
        end
        E = energies(model)
        u, i, k, at_bottom = E.fastest
        @printf("  %5d %+9.2e %+9.2e %+9.2e %+9.2e %+9.2e %7.3f (%2d,%2d)%s\n", n,
                (E.kinetic + E.potential + E.surface - total₀) / KE₀, (E.kinetic - KE₀) / KE₀,
                (E.bottom_kinetic - E₀.bottom_kinetic) / KE₀, (E.potential - E₀.potential) / KE₀,
                (E.surface - E₀.surface) / KE₀, abs(u), i, k, at_bottom ? "*" : "")
    end
    println()
end
