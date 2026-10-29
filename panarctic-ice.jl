# Where is the sea ice fast, and where is it piling up? Reads the hourly sea-ice output of a run
# (arctic_sea_ice_surface.jld2 with he = h ℵ, uᵢ, vᵢ).
#
# Usage, from the folder with the output, giving the model folder as the project:
#   julia --project=.. ../../panarctic-tools/panarctic-ice.jl
#
# Prints, every 6 hours and for the last 6 outputs: the fastest ice anywhere, the fastest ice where there is ice
# (he > 0.1 m), and the thickest ice, each with its location. Saves ice_final.png (ice volume per area and speed
# at the last output) and ice_timeseries.png.

using Oceananigans, Oceananigans.Units, Printf, CairoMakie
using Oceananigans.Grids: λnode, φnode

file = "arctic_sea_ice_surface.jld2"
he = FieldTimeSeries(file, "he")
uᵢ = FieldTimeSeries(file, "uᵢ")
vᵢ = FieldTimeSeries(file, "vᵢ")
grid = he.grid
Nx, Ny = size(grid)[1:2]
times = he.times
Nt = length(times)

ℓ = (Center(), Center(), Center())
λ = [λnode(i, j, 1, grid, ℓ...) for i in 1:Nx, j in 1:Ny]
φ = [φnode(i, j, 1, grid, ℓ...) for i in 1:Nx, j in 1:Ny]

# Ice speed at cell centres from the face velocities
function ice_speed(n)
    u = Array(interior(uᵢ[n]))[:, :, 1]
    v = Array(interior(vᵢ[n]))[:, :, 1]
    uc = [(u[i, j] + u[min(i + 1, size(u, 1)), j]) / 2 for i in 1:Nx, j in 1:Ny]
    vc = [(v[i, j] + v[i, min(j + 1, size(v, 2))]) / 2 for i in 1:Nx, j in 1:Ny]
    return replace!(sqrt.(uc .^ 2 .+ vc .^ 2), NaN => 0)
end
ice_volume(n) = replace!(Array(interior(he[n]))[:, :, 1], NaN => 0)

where(i, j) = @sprintf("(%d, %d) ≈ %.1f°E %.1f°N, %d cells from an open edge", i, j, λ[i, j], φ[i, j], min(j - 1, Ny - j))

function report(n)
    s, h = ice_speed(n), ice_volume(n)
    m₁, I₁ = findmax(s)
    m₂, I₂ = findmax(ifelse.(h .> 0.1, s, 0.0))
    m₃, I₃ = findmax(h)
    @printf("day %5.2f  fastest anywhere %.2f m/s at %s\n", times[n] / days, m₁, where(Tuple(I₁)...))
    @printf("           fastest where he > 0.1 m: %.2f m/s at %s (he there %.2f m)\n", m₂, where(Tuple(I₂)...), h[I₂])
    @printf("           thickest: he = %.2f m at %s; cells with he > 8 m: %d\n", m₃, where(Tuple(I₃)...), count(h .> 8))
end

shown = unique(vcat(1:6:Nt, max(1, Nt - 5):Nt))
foreach(report, shown)

# Time series of the extremes
fastest      = [maximum(ice_speed(n)) for n in 1:Nt]
fastest_ice  = [maximum(ifelse.(ice_volume(n) .> 0.1, ice_speed(n), 0.0)) for n in 1:Nt]
thickest     = [maximum(ice_volume(n)) for n in 1:Nt]

fig = Figure(size = (900, 600))
ax = Axis(fig[1, 1], ylabel = "ice speed [m/s]")
lines!(ax, times ./ days, fastest, label = "anywhere")
lines!(ax, times ./ days, fastest_ice, label = "where he > 0.1 m")
axislegend(ax, position = :lt)
ax = Axis(fig[2, 1], xlabel = "day", ylabel = "max he [m]")
lines!(ax, times ./ days, thickest)
save("ice_timeseries.png", fig)

# Final state: ice volume per area and speed
h, s = ice_volume(Nt), ice_speed(Nt)
s[h .< 0.01] .= NaN
_, I = findmax(h)
fig = Figure(size = (1300, 600))
ax = Axis(fig[1, 1], title = @sprintf("ice volume per area he [m], day %.2f", times[Nt] / days), aspect = DataAspect())
hm = heatmap!(ax, h, colormap = Reverse(:blues), colorrange = (0, 5), highclip = :red)
Colorbar(fig[1, 2], hm)
scatter!(ax, [I[1]], [I[2]], marker = :circle, color = :transparent, strokecolor = :red, strokewidth = 2, markersize = 20)
ax = Axis(fig[1, 3], title = "ice speed [m/s] where there is ice", aspect = DataAspect())
hm = heatmap!(ax, s, colormap = Reverse(:solar), colorrange = (0, 1))
Colorbar(fig[1, 4], hm)
save("ice_final.png", fig)
println("Saved ice_timeseries.png and ice_final.png")
