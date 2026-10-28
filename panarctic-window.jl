# A window around one (i, j): bathymetry, how far the grid axes are rotated from geographic east/north, and the
# surface flow from the daily 3-D output (arctic_ocean_full.jld2), if it is there.
#
# Usage, from the model folder:
#   julia --project=. ../panarctic-tools/panarctic-window.jl 366 31        # half-width 30 cells
#   julia --project=. ../panarctic-tools/panarctic-window.jl 366 31 20
#
# Writes window_<i>_<j>.png. The angle is that of the grid's +j axis, clockwise from geographic north:
# where it is not small, GLORYS v (northward) is not the grid's v.

using NumericalEarth, Oceananigans, Oceananigans.Units, Printf
using Oceananigans.OrthogonalSphericalShellGrids: RotatedLatitudeLongitudeGrid
using Oceananigans.Grids: λnode, φnode
using CairoMakie

i₀, j₀ = parse(Int, ARGS[1]), parse(Int, ARGS[2])
w = length(ARGS) ≥ 3 ? parse(Int, ARGS[3]) : 30

λ₁, λ₂ = 151, 212
φ₁, φ₂ = -34.8, 27
res = 0.15
Nx, Ny, Nz = round(Int, (λ₂ - λ₁) / res), round(Int, (φ₂ - φ₁) / res), 20
depth = 4000meters
z = ExponentialDiscretization(Nz, -depth, 0; scale = depth/4, mutable = true)

grid = RotatedLatitudeLongitudeGrid(CPU(); size = (Nx, Ny, Nz), longitude = (λ₁, λ₂), latitude = (φ₁, φ₂),
                                    north_pole = (145, 0), z, halo = (7, 7, 7))

bathymetry = Metadatum(:bottom_height; dataset = ETOPO2022())
h  = regrid_bathymetry(grid, bathymetry; minimum_depth = 15, interpolation_passes = 2, major_basins = 1)
hb = Array(interior(h, :, :, 1))
hb[hb .≥ 0] .= NaN

ℓ = (Center(), Center(), Center())
λ = [λnode(i, j, 1, grid, ℓ...) for i in 1:Nx, j in 1:Ny]
φ = [φnode(i, j, 1, grid, ℓ...) for i in 1:Nx, j in 1:Ny]

# Clockwise angle of the +j direction from geographic north, from neighbouring cell centres
function j_axis_angle(i, j)
    j₊, j₋ = min(j + 1, Ny), max(j - 1, 1)
    Δλ = mod(λ[i, j₊] - λ[i, j₋] + 180, 360) - 180
    east  = Δλ * cosd(φ[i, j])
    north = φ[i, j₊] - φ[i, j₋]
    return atand(east, north)
end
θ = [j_axis_angle(i, j) for i in 1:Nx, j in 1:Ny]

is = max(1, i₀ - w):min(Nx, i₀ + w)
js = max(1, j₀ - w):min(Ny, j₀ + w)

@printf("\n(i, j) = (%d, %d) ≈ %.2f°E %.2f°N, bottom %.0f m, grid +j axis %.0f° clockwise from north\n",
        i₀, j₀, λ[i₀, j₀], φ[i₀, j₀], hb[i₀, j₀], θ[i₀, j₀])
@printf("Window corners: (%d, %d) ≈ %.1f°E %.1f°N   (%d, %d) ≈ %.1f°E %.1f°N\n",
        first(is), first(js), λ[first(is), first(js)], φ[first(is), first(js)],
        last(is), last(js), λ[last(is), last(js)], φ[last(is), last(js)])

if first(js) == 1
    wet_edge = [i for i in is if isfinite(hb[i, 1])]
    if isempty(wet_edge)
        println("South open edge (j = 1) in the window: all land")
    else
        @printf("South open edge (j = 1) in the window: ocean at i = %d–%d (%d cells); grid +j axis %.0f° to %.0f° from north there\n",
                first(wet_edge), last(wet_edge), length(wet_edge), extrema(θ[wet_edge, 1])...)
    end
end

fig = Figure(size = (1500, 560))
Label(fig[0, 1:3], @sprintf("Window around (i, j) = (%d, %d), %.2f°E %.2f°N", i₀, j₀, λ[i₀, j₀], φ[i₀, j₀]), fontsize = 18)

ax = Axis(fig[1, 1][1, 1], title = "bottom height [m]", xlabel = "i", ylabel = "j", aspect = DataAspect())
hm = heatmap!(ax, is, js, hb[is, js], colormap = :deep, colorrange = (-min(500, -minimum(filter(isfinite, hb[is, js]))), 0),
              nan_color = :gray75)
Colorbar(fig[1, 1][2, 1], hm, vertical = false)
scatter!(ax, [i₀], [j₀], marker = :xcross, color = :red, markersize = 14)
contour!(ax, is, js, λ[is, js], levels = -180:5:180, color = :black, linewidth = 0.5, labels = true)
contour!(ax, is, js, φ[is, js], levels = 0:2:90, color = :black, linestyle = :dash, linewidth = 0.5, labels = true)

ax = Axis(fig[1, 2][1, 1], title = "grid +j axis, degrees clockwise from north", xlabel = "i", aspect = DataAspect())
hm = heatmap!(ax, is, js, θ[is, js], colormap = :balance, colorrange = (-90, 90))
Colorbar(fig[1, 2][2, 1], hm, vertical = false)
scatter!(ax, [i₀], [j₀], marker = :xcross, color = :black, markersize = 14)

file = "arctic_ocean_full.jld2"
if isfile(file)
    u = FieldTimeSeries(file, "u")
    v = FieldTimeSeries(file, "v")
    n = length(u.times)
    uᵢ = Array(interior(u[n], :, :, Nz))
    vᵢ = Array(interior(v[n], :, :, Nz))
    uc = [(uᵢ[i, j] + uᵢ[i+1, j]) / 2 for i in is, j in js]
    vc = [(vᵢ[i, j] + vᵢ[i, j+1]) / 2 for i in is, j in js]
    speed = sqrt.(uc .^ 2 .+ vc .^ 2)
    speed[isnan.(hb[is, js])] .= NaN

    ax = Axis(fig[1, 3][1, 1], title = @sprintf("surface speed [m/s] and direction (grid axes), day %.1f", u.times[n] / days),
              xlabel = "i", aspect = DataAspect())
    hm = heatmap!(ax, is, js, speed, colormap = Reverse(:solar), colorrange = (0, 1.5), nan_color = :gray75)
    Colorbar(fig[1, 3][2, 1], hm, vertical = false)
    step = max(1, w ÷ 10)
    pts = [Point2f(i, j) for i in is[1:step:end], j in js[1:step:end] if isfinite(hb[i, j])]
    dirs = [Vec2f(uc[i - first(is) + 1, j - first(js) + 1], vc[i - first(is) + 1, j - first(js) + 1]) for i in is[1:step:end], j in js[1:step:end] if isfinite(hb[i, j])]
    arrows2d!(ax, pts, dirs, lengthscale = 3 * step, color = :black)
    scatter!(ax, [i₀], [j₀], marker = :xcross, color = :blue, markersize = 14)

    @printf("Surface velocity at (%d, %d), day %.1f: u = %.2f, v = %.2f m/s (grid axes)\n", i₀, j₀, u.times[n] / days,
            uc[i₀ - first(is) + 1, j₀ - first(js) + 1], vc[i₀ - first(is) + 1, j₀ - first(js) + 1])
else
    println("No $file here; skipping the flow panel.")
end

path = "window_$(i₀)_$(j₀).png"
save(path, fig)
println("Saved $path")
