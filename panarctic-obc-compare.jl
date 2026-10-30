# One simple figure for the Oceananigans discussion: three boundary set-ups against GLORYS.
#
# Panels: (a) net transport into the domain through the Atlantic edge, (b) transport in the Norwegian jet band
# (i 335–352) at j = 175, (c) domain-mean free surface. Reads transport_monitor.csv from each run.
#
# Usage (any terminal, after ocean_setup_13):
#   julia --project=<model folder> panarctic-obc-compare.jl <results folder> [out.png]
# where <results folder> contains transport-1day, transport-1day-sponge and transport-1day-target.

using CairoMakie

function read_monitor(path)
    ls = readlines(path)
    header = split(ls[1], ',')
    data = [parse.(Float64, split(l, ',')) for l in ls[2:end] if !isempty(l)]
    return name -> [d[findfirst(==(name), header)] for d in data]
end

root = isempty(ARGS) ? "." : ARGS[1]
out  = length(ARGS) ≥ 2 ? ARGS[2] : joinpath(root, "obc_compare_1day.png")

runs = [("Flather (current)",        "transport-1day"),
        ("Flather + velocity sponge", "transport-1day-sponge"),
        ("Flather + target_transport", "transport-1day-target")]
cols = [:steelblue, :darkorange, :seagreen]

fig = Figure(size = (1500, 430), fontsize = 15)
axes = [Axis(fig[1, 1], title = "(a) net transport into the domain, Atlantic edge", xlabel = "hour", ylabel = "Sv"),
        Axis(fig[1, 2], title = "(b) Norwegian slope current (jet band, j = 175)", xlabel = "hour", ylabel = "Sv"),
        Axis(fig[1, 3], title = "(c) domain-mean free surface", xlabel = "hour", ylabel = "m")]
names = [("j1_model_Sv", "j1_glorys_Sv"), ("jet_model_Sv", "jet_glorys_Sv"), ("eta_model_m", "eta_glorys_m")]

for (r, (label, dir)) in enumerate(runs)
    col = read_monitor(joinpath(root, dir, "transport_monitor.csv"))
    hour = col("hour")
    for (p, (m, _)) in enumerate(names)
        lines!(axes[p], hour, col(m), color = cols[r], linewidth = 2, label = label)
    end
    if r == 1
        for (p, (_, g)) in enumerate(names)
            lines!(axes[p], hour, col(g), color = :black, linestyle = :dash, linewidth = 2, label = "GLORYS")
        end
    end
end
Legend(fig[2, 1:3], axes[1], orientation = :horizontal, framevisible = false)
save(out, fig)
println("Saved $out")
