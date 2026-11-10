# A 1° (100 km) pan-Arctic ocean nested in GLORYS with NumericalEarth's nested_ocean_model (PR #676 at 5975f67f,
# "Nest an ocean on any spherical grid"). Ocean only: no atmosphere, no sea ice. Same conformal conic grid as
# PanArcticShenanigans (λ₀ = −40°, pole inside, Atlantic edge j = 1, Pacific edge j = Ny), same dates.
#
# Set up once, in its own environment (the branch needs Oceananigans 0.117, so it must not touch panarctic-port):
#
#   cd panarctic-tools/nested
#   julia --project=. -e 'using Pkg;
#       Pkg.add(url = "https://github.com/NumericalEarth/NumericalEarth.jl", rev = "5975f67f7b3dd54519d1590dada3d7f6847ea45e");
#       Pkg.add(["Oceananigans", "CUDA", "CopernicusMarine", "JLD2"]); Pkg.precompile()'
#
# Run:  julia --project=. panarctic_nested_1deg.jl  [output_dir]
#
# What to check (printed hourly and written to nested_monitor.csv):
#   - the net transport through the Atlantic (south) and Pacific (north) edges, positive into the domain: the
#     Pacific edge should bring in about +0.7–0.9 Sv (Bering throughflow) and the Atlantic edge carry it out. A sign
#     or size far from that means the velocity rotation or the edge sampling is wrong;
#   - the domain-mean free surface: plain Flather, no volume constraint, so GLORYS's edge transports that do not
#     close show up as a drift (about 0.5 cm/day per Sv);
#   - the shelf (H < 200 m) and deep (H > 1000 m) mean η change, to compare with our 15 km runs.

using NumericalEarth
using Oceananigans
using Oceananigans.Units
using Oceananigans.Architectures: on_architecture
using Oceananigans.Operators: Δxᶜᶠᶜ, Δzᶜᶠᶜ, Azᶜᶜᶜ
using Oceananigans.OrthogonalSphericalShellGrids: ConformalConicGrid
using CopernicusMarine
using CUDA
using Dates
using Printf

output_dir = isempty(ARGS) ? joinpath(@__DIR__, "output-1deg") : ARGS[1]
mkpath(output_dir)

arch = GPU()

# ## Grid and bathymetry (as PanArcticShenanigans setup.jl, at 100 km)

spacing = 100kilometers
Nz      = 20
depth   = 5000meters
z       = ExponentialDiscretization(Nz, -depth, 0; scale = depth/4, mutable = true)

x₁, x₂ = -3280e3, 3616e3
y₁, y₂ = -4100e3, 2700e3

underlying_grid = ConformalConicGrid(arch;
                                     size = (round(Int, (x₂ - x₁) / spacing), round(Int, (y₂ - y₁) / spacing), Nz),
                                     x = (x₁, x₂),
                                     y = (y₁, y₂),
                                     central_longitude  = -40,
                                     standard_parallel  = 90,
                                     latitude_of_origin = 90,
                                     z,
                                     halo = (7, 7, 7))

Nx, Ny, _ = size(underlying_grid)
@info "Grid: $Nx × $Ny × $Nz at $(spacing / 1e3) km"
@info "Region the parent is read over: $(BoundingBox(underlying_grid))"

# ETOPO (NumericalEarth's default) rather than GEBCO: quick to regrid at 100 km
bottom_height = regrid_bathymetry(underlying_grid; minimum_depth = 15, interpolation_passes = 2, major_basins = 1)
grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom_height); active_cells_map = true)

# ## The nested ocean
#
# Only the Atlantic and Pacific edges are open (the default would open every side that touches water). Everything
# else takes the branch's defaults: oblique radiation (9 h inflow), tracer reservoirs (20 km), Flather + Chapman, a
# 2/3° sponge on T, S and the baroclinic velocity, the bathymetry blended to GLORYS's over 1/3°, and the initial
# state interpolated from GLORYS (not balanced).

dates = (DateTime(1993, 1, 1), DateTime(1993, 1, 6))

ocean = nested_ocean_simulation(grid, GLORYSDaily(); dates,
                                open_sides = (:south, :north),
                                Δt = 20minutes,
                                stop_time = 4days)

child = ocean.model.child
@info "Child grid topology: $(topology(child.grid)); open sides: $(keys(ocean.model.exchanger.bands))"

# ## Monitor: edge transports, domain-mean η, shelf and deep η

cpu_bottom = Array(interior(on_architecture(CPU(), child.grid.immersed_boundary.bottom_height)))[:, :, 1]
wet        = cpu_bottom .< 0
shelf      = wet .& (cpu_bottom .> -200)
deep       = wet .& (cpu_bottom .< -1000)

function edge_transport(v, j)
    g = v.grid
    return sum(v[i, j, k] * Δxᶜᶠᶜ(i, j, k, g) * Δzᶜᶠᶜ(i, j, k, g) for i in 1:Nx, k in 1:Nz)
end

function area_mean(η, mask, g)
    a = [Azᶜᶜᶜ(i, j, 1, g) for i in 1:Nx, j in 1:Ny]
    return sum(η[mask] .* a[mask]) / sum(a[mask])
end

η_start = Ref{Matrix{Float64}}()
monitor_file = joinpath(output_dir, "nested_monitor.csv")
open(io -> println(io, "hour,atlantic_in_Sv,pacific_in_Sv,net_in_Sv,eta_mean_cm,shelf_cm,deep_cm,max_u,max_v"), monitor_file, "w")

wall_time = Ref(time_ns())

function monitor(sim)
    v = on_architecture(CPU(), child.velocities.v)       # carries a CPU copy of the grid with the current z★ state
    η = Array(interior(on_architecture(CPU(), child.free_surface.displacement)))[:, :, 1]
    g = v.grid
    isassigned(η_start) || (η_start[] = copy(η))
    δη = η .- η_start[]

    atlantic = edge_transport(v, 1) / 1e6                 # positive northward = into the domain
    pacific  = -edge_transport(v, Ny + 1) / 1e6           # positive southward = into the domain
    η̄        = 100 * area_mean(η, wet, g)
    shelf_δη = 100 * area_mean(δη, shelf, g)
    deep_δη  = 100 * area_mean(δη, deep, g)
    umax, vmax = maximum(abs, child.velocities.u), maximum(abs, child.velocities.v)
    hour = time(sim) / hours

    @info @sprintf("%5.1f h  Atlantic in %+6.2f Sv, Pacific in %+6.2f Sv, net %+6.2f Sv; mean η %+7.2f cm; Δη shelf %+6.2f, deep %+6.2f cm; max|u|, |v| %.2f, %.2f m/s; wall %s",
                   hour, atlantic, pacific, atlantic + pacific, η̄, shelf_δη, deep_δη, umax, vmax,
                   prettytime(1e-9 * (time_ns() - wall_time[])))
    open(io -> @printf(io, "%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", hour, atlantic, pacific, atlantic + pacific,
                       η̄, shelf_δη, deep_δη, umax, vmax), monitor_file, "a")
    wall_time[] = time_ns()
    return nothing
end

add_callback!(ocean, monitor, TimeInterval(1hour))

# ## Output: surface fields and η every 3 hours

surface = (; u = child.velocities.u, v = child.velocities.v, T = child.tracers.T, S = child.tracers.S)
ocean.output_writers[:surface] = JLD2Writer(child, surface; indices = (:, :, Nz), schedule = TimeInterval(3hours),
                                            filename = joinpath(output_dir, "nested_surface"), overwrite_files = true)
ocean.output_writers[:free_surface] = JLD2Writer(child, (; η = child.free_surface.displacement),
                                                 schedule = TimeInterval(1hour),
                                                 filename = joinpath(output_dir, "nested_free_surface"),
                                                 overwrite_files = true)

run!(ocean)
