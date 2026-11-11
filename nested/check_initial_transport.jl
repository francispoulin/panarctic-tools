# Is the parent velocity rotated correctly on our grid? Build the nested ocean at 15 km (the PanArcticShenanigans
# grid) on 1993-01-01, without stepping, and print the net transport of the interpolated, rotated GLORYS velocity
# through the Atlantic and Pacific edges. PanArcticShenanigans on the same grid and day (its own rotation, in
# kernels.jl) gives, in this convention (positive into the domain):
#
#     GLORYS on the 15 km model grid:   Atlantic −0.02 Sv, Pacific +0.74 Sv
#     GLORYS on its own grid:           Atlantic −0.69 Sv, Pacific +0.91 Sv
#
# The nested model uses ETOPO and blends the seafloor to GLORYS's along the edges, so ±1 Sv is fine. Tens of Sv mean
# the rotation or the edge sampling is wrong (at 1° the first hour gave +21.7 Sv at the Atlantic edge).
#
# In the nested REPL (no GPU memory needed beyond the model; GLORYS files for 1993-01-01 are already downloaded):
#   empty!(ARGS); push!(ARGS, "15"); include("check_initial_transport.jl")       # spacing in km

include("common.jl")
using CopernicusMarine
using CUDA
using Dates

spacing_km = isempty(ARGS) ? 15 : parse(Float64, ARGS[1])
check_grid = panarctic_grid(GPU(); spacing = spacing_km * kilometers)

check_ocean = nested_ocean_simulation(check_grid, GLORYSDaily(); dates = (DateTime(1993, 1, 1), DateTime(1993, 1, 3)),
                                      open_sides = (:south, :north), Δt = 10minutes, stop_time = 1hours)

Q = edge_transports(check_ocean.model.child)
@info @sprintf("Initial state at %g km, into the domain: Atlantic %+.2f Sv, Pacific %+.2f Sv, net %+.2f Sv",
               spacing_km, Q.atlantic, Q.pacific, Q.atlantic + Q.pacific)
