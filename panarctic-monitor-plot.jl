# Plot transport_monitor.csv (written by panarctic-transport-monitor.jl): model vs GLORYS through every
# monitored row, the Norwegian jet band and the domain-mean free surface.
#
# Usage (from anywhere; the CSV path is the argument):
#   julia --project=<model folder> panarctic-monitor-plot.jl <output_dir>/transport_monitor.csv
# or in a session:  ARGS = [...]; include(...)  — or just call monitor_plot("<path>").
#
# Writes transport_monitor.png next to the CSV.

using CairoMakie

function monitor_plot(path)
    lines_ = readlines(path)
    header = split(lines_[1], ',')
    data   = [parse.(Float64, split(l, ',')) for l in lines_[2:end] if !isempty(l)]
    col(name) = [d[findfirst(==(name), header)] for d in data]
    hour = col("hour")

    rows = unique([replace(h, r"_(model|glorys)_Sv$" => "") for h in header if occursin(r"^j\d+_", h)])
    panels = vcat(rows, ["jet"])

    fig = Figure(size = (1200, 300 * cld(length(panels) + 1, 2)))
    for (p, name) in enumerate(panels)
        ax = Axis(fig[cld(p, 2), mod1(p, 2)], title = name == "jet" ? "jet band (i 335–352) at j = 175" :
                  "row $(name[2:end]) (+ = into domain at edges, towards +j inside)", ylabel = "Sv")
        lines!(ax, hour, col(name * "_model_Sv"), label = "model")
        lines!(ax, hour, col(name * "_glorys_Sv"), label = "GLORYS", linestyle = :dash)
        p == 1 && axislegend(ax, position = :lt)
    end
    p = length(panels) + 1
    ax = Axis(fig[cld(p, 2), mod1(p, 2)], title = "domain-mean η", ylabel = "m", xlabel = "hour")
    lines!(ax, hour, col("eta_model_m"), label = "model")
    lines!(ax, hour, col("eta_glorys_m"), label = "GLORYS", linestyle = :dash)

    out = joinpath(dirname(path), "transport_monitor.png")
    save(out, fig)
    println("Saved $out")
    return fig
end

isempty(ARGS) || monitor_plot(ARGS[1])
