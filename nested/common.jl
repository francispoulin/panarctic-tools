# Shared by the nested pan-Arctic scripts: the PanArcticShenanigans grid at any spacing, and an hourly monitor of the
# edge transports and of η (domain mean, shelf, deep).

using NumericalEarth
using Oceananigans
using Oceananigans.Units
using Oceananigans.Architectures: on_architecture
using Oceananigans.Operators: Δxᶜᶠᶜ, Δzᶜᶠᶜ, Azᶜᶜᶜ
using Oceananigans.OrthogonalSphericalShellGrids: ConformalConicGrid
using Printf

# The PanArcticShenanigans conformal conic grid (λ₀ = −40°, pole inside; Atlantic edge j = 1, Pacific edge j = Ny),
# with ETOPO bathymetry (NumericalEarth's default), minimum depth 15 m and only the largest basin kept
function panarctic_grid(arch; spacing, Nz = 20, depth = 5000meters)
    z = ExponentialDiscretization(Nz, -depth, 0; scale = depth/4, mutable = true)
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
    @info "Grid: $Nx × $Ny × $Nz at $(spacing / 1e3) km; parent region $(BoundingBox(underlying_grid))"

    bottom_height = regrid_bathymetry(underlying_grid; minimum_depth = 15, interpolation_passes = 2, major_basins = 1)
    return ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom_height); active_cells_map = true)
end

# Net transport through the Atlantic (south, j = 1) and Pacific (north, j = Ny + 1) edges, positive into the domain (Sv).
# PanArcticShenanigans prints "GLORYS at the start, model grid: Atlantic a Sv in, Pacific p Sv out": in this convention
# those are (a, −p).
function edge_transports(child)
    Nx, Ny, Nz = size(child.grid)
    v = on_architecture(CPU(), child.velocities.v)       # with a CPU copy of the grid in its current z★ state
    Q(j) = sum(v[i, j, k] * Δxᶜᶠᶜ(i, j, k, v.grid) * Δzᶜᶠᶜ(i, j, k, v.grid) for i in 1:Nx, k in 1:Nz)
    return (atlantic = Q(1) / 1e6, pacific = -Q(Ny + 1) / 1e6)
end

# Hourly: net transport through the Atlantic (south) and Pacific (north) edges, positive into the domain; domain-mean
# η; the change in η since the first call averaged over the shelf (H < 200 m) and the deep ocean (H > 1000 m).
# Written to <output_dir>/nested_monitor.csv.
function add_nested_monitor!(simulation, child, output_dir; schedule = TimeInterval(1hours))
    Nx, Ny, Nz = size(child.grid)
    bottom = Array(interior(bottom_height_field(child.grid)))[:, :, 1]   # in 0.117 the immersed boundary holds a bare array
    wet    = bottom .< 0
    shelf  = wet .& (bottom .> -200)
    deep   = wet .& (bottom .< -1000)

    function area_mean(a, mask, g)
        A = [Azᶜᶜᶜ(i, j, 1, g) for i in 1:Nx, j in 1:Ny]
        return sum(a[mask] .* A[mask]) / sum(A[mask])
    end

    file = joinpath(output_dir, "nested_monitor.csv")
    open(io -> println(io, "hour,atlantic_in_Sv,pacific_in_Sv,net_in_Sv,eta_mean_cm,shelf_change_cm,deep_change_cm,max_u,max_v"),
         file, "w")

    η₀   = Ref{Matrix{Float64}}()
    wall = Ref(time_ns())

    function nested_monitor(sim)
        v = on_architecture(CPU(), child.velocities.v)     # for the grid metrics in the current z★ state
        η = Array(interior(on_architecture(CPU(), child.free_surface.displacement)))[:, :, 1]
        isassigned(η₀) || (η₀[] = copy(η))
        δη = η .- η₀[]

        atlantic, pacific = edge_transports(child)
        η̄        = 100 * area_mean(η, wet, v.grid)
        shelf_δη = 100 * area_mean(δη, shelf, v.grid)
        deep_δη  = 100 * area_mean(δη, deep, v.grid)
        umax, vmax = maximum(abs, child.velocities.u), maximum(abs, child.velocities.v)
        model_hour = time(sim) / hours

        @info @sprintf("%6.1f h  Atlantic in %+6.2f Sv, Pacific in %+6.2f Sv, net %+6.2f Sv; mean η %+7.2f cm; Δη shelf %+6.2f, deep %+6.2f cm; max|u|, |v| %.2f, %.2f m/s; wall %s",
                       model_hour, atlantic, pacific, atlantic + pacific, η̄, shelf_δη, deep_δη, umax, vmax,
                       prettytime(1e-9 * (time_ns() - wall[])))
        open(io -> @printf(io, "%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", model_hour, atlantic, pacific,
                           atlantic + pacific, η̄, shelf_δη, deep_δη, umax, vmax), file, "a")
        wall[] = time_ns()
        return nothing
    end

    add_callback!(simulation, nested_monitor, schedule; name = :nested_monitor)
    return file
end

# Surface fields every 3 hours and η every hour, from the child ocean
function add_nested_output!(simulation_or_ocean, child, output_dir)
    Nz = size(child.grid, 3)
    writers = simulation_or_ocean.output_writers
    surface = (; u = child.velocities.u, v = child.velocities.v, T = child.tracers.T, S = child.tracers.S)
    writers[:surface] = JLD2Writer(child, surface; indices = (:, :, Nz), schedule = TimeInterval(3hours),
                                   filename = joinpath(output_dir, "nested_surface"), overwrite_files = true)
    writers[:free_surface] = JLD2Writer(child, (; η = child.free_surface.displacement), schedule = TimeInterval(1hours),
                                        filename = joinpath(output_dir, "nested_free_surface"), overwrite_files = true)
    return nothing
end

# Per-column transport through each edge (m³/s, positive into the domain) and the column depth there, with the
# column's longitude and latitude, written to <file> as CSV. Compare with port_edge_columns.jl (the same columns from
# PanArcticShenanigans) using compare_edge_columns.jl.
function write_edge_columns(child, file)
    Nx, Ny, Nz = size(child.grid)
    v = on_architecture(CPU(), child.velocities.v)
    g = v.grid
    open(file, "w") do io
        println(io, "edge,i,lon,lat,depth_m,transport_m3s")
        for (edge, j, sign) in (("atlantic", 1, 1), ("pacific", Ny + 1, -1)), i in 1:Nx
            U = sign * sum(v[i, j, k] * Δxᶜᶠᶜ(i, j, k, g) * Δzᶜᶠᶜ(i, j, k, g) for k in 1:Nz)
            H = max(Oceananigans.Grids.static_column_depthᶜᶠᵃ(i, j, g), 0.0)
            λ = Oceananigans.Grids.λnode(i, j, 1, g, Center(), Face(), Center())
            φ = Oceananigans.Grids.φnode(i, j, 1, g, Center(), Face(), Center())
            @printf(io, "%s,%d,%.4f,%.4f,%.2f,%.6e\n", edge, i, λ, φ, H, U)
        end
    end
    @info "Wrote $file"
    return file
end
