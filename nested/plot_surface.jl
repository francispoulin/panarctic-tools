# Maps of a nested run's surface output, first and last snapshot side by side, and the domain-mean η over time.
# Works for the branch's Arctic example (arctic_surface.jld2: T, S, u, v, η, h, ℵ every 2 days) and for our scripts
# (nested_surface.jld2: T, S, u, v every 3 hours; η from nested_free_surface.jld2).
#
#   julia --project=. plot_surface.jl output-arctic-example/arctic_surface.jld2
#
# Needs CairoMakie in the nested environment (once: julia --project=. -e 'using Pkg; Pkg.add("CairoMakie")').
# Writes <file>_maps.png and <file>_eta_mean.png next to the file. Maps are in grid index space (i, j).

using Oceananigans
using JLD2
using CairoMakie
using Printf

file = ARGS[1]
stem = splitext(file)[1]

names = filter(n -> n in ("T", "S", "u", "v", "η", "h", "ℵ"), jldopen(f -> collect(keys(f["timeseries"])), file))
series = Dict(n => FieldTimeSeries(file, n; backend = InMemory()) for n in names)
times = series[first(names)].times
Nt = length(times)
@info "$(basename(file)): $(join(names, ", ")); $Nt snapshots from $(prettytime(times[1])) to $(prettytime(times[end]))"

surface(fts, n) = (a = interior(fts[n]); size(a, 3) > 1 ? a[:, :, end] : a[:, :, 1])
masked(a) = map(x -> x == 0 ? NaN : x, a)          # land (exactly zero in the output) shown blank

ranges = Dict("T" => (-2, 12), "S" => (30, 35.5), "u" => (-0.3, 0.3), "v" => (-0.3, 0.3), "η" => (-1, 1),
              "h" => (0, 4), "ℵ" => (0, 1))
colormaps = Dict("T" => :thermal, "S" => :haline, "u" => :balance, "v" => :balance, "η" => :balance, "h" => :ice,
                 "ℵ" => :ice)

fig = Figure(size = (900, 330 * length(names)))
for (row, n) in enumerate(names)
    for (col, k) in enumerate((1, Nt))
        a = surface(series[n], k)
        ax = Axis(fig[row, col]; title = @sprintf("%s at %s", n, prettytime(times[k])), aspect = DataAspect())
        hm = heatmap!(ax, n in ("h", "ℵ") ? a : masked(a); colormap = colormaps[n], colorrange = ranges[n])
        col == 2 && Colorbar(fig[row, 3], hm)
    end
end
save(stem * "_maps.png", fig)
@info "Wrote $(stem)_maps.png"

if "η" in names
    η̄ = [sum(filter(!isnan, masked(surface(series["η"], k)))) / count(!isnan, masked(surface(series["η"], k)))
          for k in 1:Nt]
    fig = Figure(size = (700, 300))
    ax = Axis(fig[1, 1]; xlabel = "days", ylabel = "cm", title = "mean η over the wet cells (unweighted)")
    lines!(ax, times ./ 86400, 100 .* η̄)
    save(stem * "_eta_mean.png", fig)
    @info @sprintf("Wrote %s_eta_mean.png; mean η %+.2f cm → %+.2f cm", stem, 100η̄[1], 100η̄[end])
end
