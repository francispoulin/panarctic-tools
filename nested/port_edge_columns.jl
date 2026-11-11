# The per-column GLORYS transport through the Atlantic and Pacific edges as PanArcticShenanigans sees it (GLORYS
# rotated and interpolated to the model's faces, kernels.jl), on the start day. Include it in a panarctic-port REPL
# in which panarctic.jl has run at 15 km (any config with 1993-01-01; the model can still be loaded or the run
# finished):
#
#   include("../panarctic-tools/nested/port_edge_columns.jl")
#
# Writes port_edge_columns.csv in the current folder, in the same format as write_edge_columns in common.jl.

using Printf
using Oceananigans.Grids: static_column_depthᶜᶠᵃ, λnode, φnode

let file = "port_edge_columns.csv", t = Oceananigans.Units.Time(vᶠ.times[start_index])
    CUDA.@allowscalar open(file, "w") do io
        println(io, "edge,i,lon,lat,depth_m,transport_m3s")
        for (edge, e, sign) in (("atlantic", 1, 1), ("pacific", 2, -1)), i in 1:Nx
            j = edge_rows[e]
            U = sign * vertical_integral(i, j, grid, vᶠ, t, Δzᶜᶠᶜ, Center(), Face(), Center()) * Δx_edge[i, e]
            @printf(io, "%s,%d,%.4f,%.4f,%.2f,%.6e\n", edge, i, λnode(i, j, 1, grid, Center(), Face(), Center()),
                    φnode(i, j, 1, grid, Center(), Face(), Center()), H_edge[i, e], U)
        end
    end
    @info "Wrote $file; GLORYS on the model grid at day start_index: Atlantic $(round(atlantic_glorys[start_index] / 1e6, digits = 2)) Sv in, Pacific $(round(-pacific_glorys[start_index] / 1e6, digits = 2)) Sv in"
end
