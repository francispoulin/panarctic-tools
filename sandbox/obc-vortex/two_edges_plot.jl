# Plot the two-edge runs (two_edges.jl):  julia --project=$OCEAN_HPC_ROOT/viz two_edges_plot.jl
using CairoMakie, Serialization
cases   = (:base, :datum_both, :datum_south, :parent_mismatch)
variant = isempty(ARGS) ? "flat" : ARGS[1]
prefix  = variant == "flat" ? "two_edges_" : "two_edges_$(variant)_"
methods = variant == "shallow" ? (:flather, :total_zero, :hybrid, :hybrid_H) : (:flather, :edge_net, :total_zero, :hybrid)
colours = Dict(:flather => :steelblue, :edge_net => :seagreen, :total_zero => :darkorange, :hybrid => :black,
               :hybrid_H => :purple)
labels  = Dict(:flather => "Flather", :edge_net => "per-edge net (M5b)", :total_zero => "total zero (NEMO/ROMS)",
               :hybrid => "hybrid (η shift, ∝ √(gH))", :hybrid_H => "hybrid (uniform velocity, ∝ H)")
fig = Figure(size = (1500, 1250), fontsize = 14)
for (row, c) in enumerate(cases)
    axs = [Axis(fig[row, 1], ylabel = "Sv", title = "$c: outflow, south edge"),
           Axis(fig[row, 2], ylabel = "Sv", title = "$c: inflow, strait"),
           Axis(fig[row, 3], ylabel = "cm", title = "$c: change in basin-mean η")]
    for m in methods
        r = deserialize("results/$(prefix)$(m)_$(c).jls")
        d = r.times ./ 86400
        lines!(axs[1], d, -r.south ./ 1e6, color = colours[m], label = labels[m])
        lines!(axs[2], d, -r.north ./ 1e6, color = colours[m])
        lines!(axs[3], d, 100 .* (r.mean_η .- r.mean_η[1]), color = colours[m])
    end
    row == length(cases) && foreach(a -> a.xlabel = "day", axs)
end
Legend(fig[length(cases) + 1, 1:3], content(fig[1, 1]), orientation = :horizontal, framevisible = false)
out = variant == "flat" ? "figures/two_edges.png" : "figures/two_edges_$(variant).png"
mkpath("figures"); save(out, fig)
println("saved $out")
