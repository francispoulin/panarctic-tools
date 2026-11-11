# Where along each edge do the nested model's and PanArcticShenanigans' initial transports differ?
#
#   julia --project=. compare_edge_columns.jl nested_edge_columns.csv port_edge_columns.csv
#
# Prints the net per edge, and the cumulative difference along each edge in eighths, so that a difference that grows
# where the grid is most rotated (the ends of each edge) points to the rotation, one that sits where the depths differ
# points to the bathymetry, and one spread evenly points to sampling. Writes edge_columns.png (cumulative transport
# and depth along each edge, both models).

using CairoMakie
using Printf

function read_columns(file)
    rows = [split(l, ",") for l in readlines(file)[2:end]]
    return Dict(edge => (i     = [parse(Int, r[2]) for r in rows if r[1] == edge],
                         lon   = [parse(Float64, r[3]) for r in rows if r[1] == edge],
                         lat   = [parse(Float64, r[4]) for r in rows if r[1] == edge],
                         depth = [parse(Float64, r[5]) for r in rows if r[1] == edge],
                         U     = [parse(Float64, r[6]) for r in rows if r[1] == edge]) for edge in ("atlantic", "pacific"))
end

nested, port = read_columns(ARGS[1]), read_columns(ARGS[2])

fig = Figure(size = (1100, 650))
for (col, edge) in enumerate(("atlantic", "pacific"))
    n, p = nested[edge], port[edge]
    length(n.i) == length(p.i) || error("$edge: the two files have different numbers of columns")
    @printf("\n%s edge: net nested %+.2f Sv, port %+.2f Sv, difference %+.2f Sv; gross nested %.1f Sv, port %.1f Sv\n",
            edge, sum(n.U) / 1e6, sum(p.U) / 1e6, (sum(n.U) - sum(p.U)) / 1e6, sum(abs, n.U) / 2e6, sum(abs, p.U) / 2e6)
    @printf("   columns where only one model is wet: %d; mean |depth difference| where both are wet: %.0f m\n",
            count(xor.(n.depth .> 0, p.depth .> 0)), sum(abs.(n.depth .- p.depth)[(n.depth .> 0) .& (p.depth .> 0)]) /
            max(1, count((n.depth .> 0) .& (p.depth .> 0))))
    N = length(n.i)
    for s in 1:8
        r = (1 + (s - 1) * N ÷ 8):(s == 8 ? N : s * N ÷ 8)
        @printf("   i %3d–%3d (%7.1f°E %5.1f°N to %7.1f°E %5.1f°N): nested %+6.2f, port %+6.2f, difference %+6.2f Sv\n",
                first(r), last(r), n.lon[first(r)], n.lat[first(r)], n.lon[last(r)], n.lat[last(r)],
                sum(n.U[r]) / 1e6, sum(p.U[r]) / 1e6, (sum(n.U[r]) - sum(p.U[r])) / 1e6)
    end
    ax1 = Axis(fig[1, col]; title = "$edge edge: cumulative transport into the domain", ylabel = "Sv")
    lines!(ax1, n.i, cumsum(n.U) ./ 1e6; label = "nested (PR 676)")
    lines!(ax1, p.i, cumsum(p.U) ./ 1e6; label = "PanArcticShenanigans")
    axislegend(ax1; position = :lt)
    ax2 = Axis(fig[2, col]; xlabel = "i", ylabel = "column depth (m)")
    lines!(ax2, n.i, n.depth)
    lines!(ax2, p.i, p.depth)
end
save("edge_columns.png", fig)
@info "Wrote edge_columns.png"
