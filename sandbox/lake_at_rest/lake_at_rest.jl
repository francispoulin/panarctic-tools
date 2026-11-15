# Lake at rest over a cliff: horizontally uniform T(z), S(z), no flow, no forcing, uniform η. The exact solution is
# no motion, so any horizontal pressure-gradient force and any velocity that grows is a numerical error.
#
# Geometry and vertical grid as the pan-Arctic model at 15 km: one cliff from 2837 m to 1863 m across one cell
# (the (306, 170) point of the 2005 run, where a 3 m/s jet sat in the bottom partial cell at z ≈ −1988 m),
# 20 exponential levels to 5000 m, z-star, TEOS-10 with ρ₀ = 1020, f-plane at 72°N.
#
# Cases: partial vs whole bottom cells, η = 0 vs η = −0.47 m (the 2005 run's mean η).
#
# Part 1, at t = 0: the force −∂x p′ at every active u face, split into
#   along the levels:  −δx p′ / Δx
#   z-star chain rule: +(∂x z)(∂z p′)          (zero when η = 0: the levels are flat)
# and again with the immersed cells below each column's bottom given that column's bottom T and S instead of the
# T = S = 0 the model holds there (does the bottom cell's force come from the immersed values?).
# Part 2: two days of stepping; the largest |u| and where it is.
#
#   julia --project=<PanArcticShenanigans>/nested lake_at_rest.jl [profile.csv] [Nz]

using NumericalEarth, Oceananigans, Oceananigans.Units, Printf
using Oceananigans.Operators: ∂xᶠᶜᶜ, δxᶠᶜᶜ, Δx⁻¹ᶠᶜᶜ
using Oceananigans.Grids: peripheral_node, inactive_node, znode
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Models.NonhydrostaticModels: update_hydrostatic_pressure!

TEOS10EquationOfState = NumericalEarth.Oceans.TEOS10EquationOfState

profile_file = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(@__DIR__, "glorys_profile_0E_72N.csv")
Nz           = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 20

# ## The profile, as functions of z (linear between GLORYS levels, constant above and below)

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

# ## Grid and model

Nx, Ny  = 16, 3
spacing = 15kilometers
depth   = 5000
deep_floor, shallow_floor = -2837, -1863
cliff_x = Nx / 2 * spacing
bottom(x, y) = x < cliff_x ? deep_floor : shallow_floor

function lake(bottom_kind, η₀)
    z = ExponentialDiscretization(Nz, -depth, 0; scale = depth / 4, mutable = true)
    underlying_grid = RectilinearGrid(CPU(); size = (Nx, Ny, Nz), x = (0, Nx * spacing), y = (0, Ny * spacing), z,
                                      halo = (7, 7, 7), topology = (Bounded, Bounded, Bounded))
    Bottom = bottom_kind == :partial ? PartialCellBottom : GridFittedBottom
    grid = ImmersedBoundaryGrid(underlying_grid, Bottom(bottom))

    model = HydrostaticFreeSurfaceModel(grid;
                                        momentum_advection = WENOVectorInvariant(),
                                        tracer_advection = WENO(order = 5),
                                        buoyancy = SeawaterBuoyancy(equation_of_state = TEOS10EquationOfState(reference_density = 1020)),
                                        coriolis = FPlane(latitude = 72),
                                        free_surface = SplitExplicitFreeSurface(grid; substeps = 30),
                                        tracers = (:T, :S))

    set!(model; η = η₀)              # first, so that z-star has moved the levels before T and S are set at them
    set!(model; T = Tᵢ, S = Sᵢ)
    return model
end

# ## Part 1: the pressure-gradient force at t = 0

bottom_cell(i, j, k, grid) = inactive_node(i - 1, j, k - 1, grid, Center(), Center(), Center()) |
                             inactive_node(i,     j, k - 1, grid, Center(), Center(), Center())

function force_report(model, label)
    grid, p = model.grid, model.pressure.pHY′
    largest = (total = (0.0, 0, 0), along = (0.0, 0, 0), chain = (0.0, 0, 0))
    largest_interior = 0.0
    for i in 2:Nx, j in 1:Ny, k in 1:Nz
        peripheral_node(i, j, k, grid, Face(), Center(), Center()) && continue
        total = -∂xᶠᶜᶜ(i, j, k, grid, p)
        along = -δxᶠᶜᶜ(i, j, k, grid, p) * Δx⁻¹ᶠᶜᶜ(i, j, k, grid)
        chain = total - along
        for (name, value) in pairs((; total, along, chain))
            abs(value) > abs(largest[name][1]) && (largest = merge(largest, NamedTuple{(name,)}(((value, i, k),))))
        end
        bottom_cell(i, j, k, grid) || (largest_interior = max(largest_interior, abs(total)))
    end
    a, i, k = largest.total
    @printf("%-38s total %9.2e m/s² at (i, k) = (%2d, %2d), z = %7.1f m%s → %6.3f m/s after a day\n", label, a, i, k,
            znode(i, 1, k, grid, Face(), Center(), Center()), bottom_cell(i, 1, k, grid) ? " (bottom cell)" : "", abs(a) * day)
    @printf("%-38s along the levels %9.2e (k = %2d), chain rule %9.2e (k = %2d), largest away from the bottom %9.2e\n", "",
            largest.along[1], largest.along[3], largest.chain[1], largest.chain[3], largest_interior)
end

function fill_immersed_tracers!(model)
    grid = model.grid
    for c in (model.tracers.T, model.tracers.S)
        data = interior(c)
        for i in 1:Nx, j in 1:Ny
            kᵇ = findfirst(k -> !inactive_node(i, j, k, grid, Center(), Center(), Center()), 1:Nz)
            isnothing(kᵇ) || (data[i, j, 1:kᵇ-1] .= data[i, j, kᵇ])
        end
        fill_halo_regions!(c)
    end
    update_hydrostatic_pressure!(model.pressure.pHY′, CPU(), grid, model.buoyancy, model.tracers)
    fill_halo_regions!(model.pressure.pHY′; only_local_halos = true)
    return nothing
end

cases = [(bottom_kind, η₀) for bottom_kind in (:partial, :whole) for η₀ in (0.0, -0.47)]

println("\n## Part 1: pressure-gradient force at t = 0 (exact: zero), Nz = $Nz")
for (bottom_kind, η₀) in cases
    model = lake(bottom_kind, η₀)
    force_report(model, @sprintf("%s cells, η = %+.2f m", bottom_kind, η₀))
    fill_immersed_tracers!(model)
    force_report(model, "  same, immersed cells filled")
end

# ## Part 2: two days of stepping

function largest_speed(model)
    grid = model.grid
    u = interior(model.velocities.u)
    best = (0.0, 0, 0)
    for i in 2:Nx, j in 1:Ny, k in 1:Nz
        peripheral_node(i, j, k, grid, Face(), Center(), Center()) && continue
        abs(u[i, j, k]) > abs(best[1]) && (best = (u[i, j, k], i, k))
    end
    return best
end

println("\n## Part 2: largest |u| during two days (exact: zero)")
for (bottom_kind, η₀) in cases
    model = lake(bottom_kind, η₀)
    line = @sprintf("%-26s", @sprintf("%s cells, η = %+.2f m", bottom_kind, η₀))
    for n in 1:2*24*6
        time_step!(model, 10minutes)
        if n % 36 == 0
            u, i, k = largest_speed(model)
            line *= @sprintf("  %3d h: %7.4f (%2d,%2d)", n ÷ 6, abs(u), i, k)
        end
    end
    println(line)
end
println("speeds in m/s, (i, k) of the u face; the cliff is the face i = $(Nx ÷ 2 + 1)")
