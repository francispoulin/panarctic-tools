# Model current speed against GLORYS around one (i, j), on one day, at the level of the model's surface output.
# Reads the grid from the run's own output, so it works for any grid the run used.
#
# Usage, from a run's output folder, giving the model folder as the project:
#   julia --project=<model folder> <tools folder>/panarctic-speed-vs-glorys.jl 274 131 3
#   julia --project=<model folder> <tools folder>/panarctic-speed-vs-glorys.jl 274 131 3 40   # half-width 40 cells
#
# Uses arctic_ocean_surface.jld2 (speed 𝒱 at level Nz - 2) and, if present, arctic_sea_ice_surface.jld2.
# GLORYS is the daily mean whose centre (noon) is nearest the chosen model time; the model value is the
# snapshot at that time. Writes speed_<i>_<j>_day<d>.png.

using NumericalEarth, Oceananigans, Oceananigans.Units, Printf, Dates, CairoMakie
using CopernicusMarine
using Oceananigans.Grids: λnode, φnode, znode

i₀, j₀ = parse(Int, ARGS[1]), parse(Int, ARGS[2])
day    = parse(Float64, ARGS[3])
w      = length(ARGS) ≥ 4 ? parse(Int, ARGS[4]) : 30
start  = DateTime(1993, 1, 1)

𝒱 = FieldTimeSeries("arctic_ocean_surface.jld2", "𝒱")
grid = 𝒱.grid
Nx, Ny, Nz = size(grid)
k = Nz - 2                                    # the level the run's surface output is written at
n = argmin(abs.(𝒱.times ./ days .- day))
t = 𝒱.times[n] / days

is = max(1, i₀ - w):min(Nx, i₀ + w)
js = max(1, j₀ - w):min(Ny, j₀ + w)
ℓ  = (Center(), Center(), Center())

model_speed = Array(interior(𝒱[n]))[is, js, 1]

date   = start + Day(floor(Int, t - 1/2))
region = BoundingBox(longitude = (-180, 180), latitude = (43, 90))
inpainting = NumericalEarth.DataWrangling.NearestNeighborInpainting(100)
uᴳ = Field(Metadatum(:u_velocity; date, dataset = GLORYSDaily(), region), grid; inpainting)
vᴳ = Field(Metadatum(:v_velocity; date, dataset = GLORYSDaily(), region), grid; inpainting)
glorys_speed = sqrt.(Array(interior(uᴳ, is, js, k)) .^ 2 .+ Array(interior(vᴳ, is, js, k)) .^ 2)

land = model_speed .== 0
glorys_speed[land] .= NaN
model_speed[land]  .= NaN

ice = nothing
if isfile("arctic_sea_ice_surface.jld2")
    he = FieldTimeSeries("arctic_sea_ice_surface.jld2", "he")
    nᵢ = argmin(abs.(he.times ./ days .- t))
    ice = Array(interior(he[nᵢ]))[is, js, 1]
end

pick(a) = a[i₀ - first(is) + 1, j₀ - first(js) + 1]
@printf("\n(i, j) = (%d, %d) ≈ %.2f°E %.2f°N, level k = %d (z ≈ %.0f m)\n", i₀, j₀,
        λnode(i₀, j₀, k, grid, ℓ...), φnode(i₀, j₀, k, grid, ℓ...), k, znode(i₀, j₀, k, grid, ℓ...))
@printf("Model snapshot at day %.2f; GLORYS daily mean of %s\n", t, Dates.format(date, "d u yyyy"))
@printf("Speed at the point:     model %.2f m/s, GLORYS %.2f m/s\n", pick(model_speed), pick(glorys_speed))
@printf("Largest in the window:  model %.2f m/s, GLORYS %.2f m/s\n",
        maximum(filter(isfinite, model_speed)), maximum(filter(isfinite, glorys_speed)))
isnothing(ice) || @printf("Ice volume per area at the point: %.2f m\n", pick(ice))

fig = Figure(size = (1300, 600))
Label(fig[0, 1:2], @sprintf("Current speed at z ≈ %.0f m around (%d, %d), %.1f°E %.1f°N", znode(i₀, j₀, k, grid, ℓ...),
                            i₀, j₀, λnode(i₀, j₀, k, grid, ℓ...), φnode(i₀, j₀, k, grid, ℓ...)), fontsize = 18)
for (col, (title, a)) in enumerate((("model, day " * @sprintf("%.2f", t), model_speed),
                                    ("GLORYS, mean of " * Dates.format(date, "d u yyyy"), glorys_speed)))
    ax = Axis(fig[1, col][1, 1], title = title, xlabel = "i", ylabel = "j", aspect = DataAspect())
    hm = heatmap!(ax, is, js, a, colormap = Reverse(:solar), colorrange = (0, 1.5), nan_color = :gray75)
    Colorbar(fig[1, col][2, 1], hm, vertical = false, label = "speed [m/s]")
    isnothing(ice) || contour!(ax, is, js, ice, levels = [0.1], color = :cyan, linewidth = 1.5)
    scatter!(ax, [i₀], [j₀], marker = :xcross, color = :blue, markersize = 14)
end
path = @sprintf("speed_%d_%d_day%.0f.png", i₀, j₀, t)
save(path, fig)
println("Saved $path (cyan line: ice edge, he = 0.1 m)")
