# Vertical sections through the Baffin Bay hot spot from arctic_ocean_full.jld2.
# Run from the repository folder:  julia --project=. ~/Software/panarctic-section.jl
# Writes panarctic_section.png.

using Oceananigans, CairoMakie, Printf
using Oceananigans.ImmersedBoundaries: immersed_inactive_node

file = "arctic_ocean_full.jld2"
T = FieldTimeSeries(file, "T")
S = FieldTimeSeries(file, "S")
u = FieldTimeSeries(file, "u")
v = FieldTimeSeries(file, "v")

grid = Oceananigans.Architectures.on_architecture(CPU(), T.grid)
Nx, Ny, Nz = size(grid)

# Cell-centre depths of the 20 levels (from the probe; same grid as the script)
zc = [-3630.7, -2959.0, -2409.1, -1958.9, -1590.3, -1288.5, -1041.4, -839.1, -673.5, -537.9,
      -426.8, -335.9, -261.5, -200.6, -150.7, -109.8, -76.4, -49.0, -26.6, -8.3]

wet = [!immersed_inactive_node(i, j, k, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny, k in 1:Nz]
masked(a) = (b = copy(a); b[.!wet] .= NaN; b)

# Average a face field to cell centres along `dim` (if it has one more point than the centres)
function to_centres(a, dim, N)
    size(a, dim) == N && return a
    lo = selectdim(a, dim, 1:N)
    hi = selectdim(a, dim, 2:N+1)
    return (lo .+ hi) ./ 2
end

function speed(n)
    uc = to_centres(Array(interior(u[n])), 1, Nx)
    vc = to_centres(Array(interior(v[n])), 2, Ny)
    return masked(sqrt.(uc .^ 2 .+ vc .^ 2))
end

tracer(fts, n) = masked(Array(interior(fts[n])))

times = T.times
ns = unique([1, min(2, length(times)), length(times)])   # t = 0, 1 h, last
hot, ref = (134, 136), (140, 155)
j₀, is = hot[2], 104:164                                  # zonal section through the hot spot

fig = Figure(size = (1500, 300 * length(ns) + 450))
panels = ((:T, "T [°C]", :thermal, (-2, 4)), (:S, "S [psu]", :haline, (32, 35)), (:speed, "speed [m/s]", :solar, (0, 3)))

for (r, n) in enumerate(ns)
    fields = (T = tracer(T, n), S = tracer(S, n), speed = speed(n))
    for (c, (name, label, cmap, crange)) in enumerate(panels)
        ax = Axis(fig[r, 2c-1]; title = @sprintf("%s, j = %d, t = %.0f h", label, j₀, times[n] / 3600),
                  xlabel = "i", ylabel = "z [m]")
        hm = heatmap!(ax, collect(is), zc, fields[name][is, j₀, :]; colormap = cmap, colorrange = crange, nan_color = :gray80)
        vlines!(ax, [hot[1]], color = :black, linestyle = :dash)
        ylims!(ax, -2600, 0)
        Colorbar(fig[r, 2c], hm)
    end
end

# Profiles at the hot spot and the reference point, first and last output
r = length(ns) + 1
axT = Axis(fig[r, 1:2]; title = "T profiles", xlabel = "T [°C]", ylabel = "z [m]")
axS = Axis(fig[r, 3:4]; title = "S profiles", xlabel = "S [psu]", ylabel = "z [m]")
for (point, colour) in ((hot, :red), (ref, :blue)), (n, style) in ((ns[1], :solid), (ns[end], :dash))
    i, j = point
    lab = @sprintf("(%d, %d), t = %.0f h", i, j, times[n] / 3600)
    lines!(axT, tracer(T, n)[i, j, :], zc; color = colour, linestyle = style, label = lab)
    lines!(axS, tracer(S, n)[i, j, :], zc; color = colour, linestyle = style, label = lab)
end
ylims!(axT, -2600, 0); ylims!(axS, -2600, 0)
axislegend(axS; position = :lb)

save("panarctic_section.png", fig)
println("saved panarctic_section.png")
