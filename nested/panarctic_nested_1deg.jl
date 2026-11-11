# Stage 1: a 1° (100 km) pan-Arctic ocean nested in GLORYS with NumericalEarth's nested_ocean_model (PR #676,
# branch ss/nested-ocean-model). Ocean only: no atmosphere, no sea ice. The PanArcticShenanigans grid (common.jl), the
# same dates as our runs.
#
# Set up once, in its own environment (the branch needs Oceananigans 0.117, so it must not touch panarctic-port):
#
#   cd panarctic-tools/nested
#   julia --project=. -e 'using Pkg; Pkg.Registry.update();
#       Pkg.add(url = "https://github.com/NumericalEarth/NumericalEarth.jl", rev = "768ed1d34d1804bc26a927de5b5fa56e258a7e57");
#       Pkg.add(["Oceananigans", "CUDA", "CopernicusMarine", "JLD2"]); Pkg.precompile()'
#
# (768ed1d3 is the head of the PR on 2026-11-10 and includes 2876238c, the fix for the bathymetry blend on curvilinear
# grids that 5975f67f lacks: "objects of type Nothing are not callable" in interpolate_bottom_height!.)
#
# Run:  julia --project=. panarctic_nested_1deg.jl  [output_dir]
#
# What to check (printed hourly and written to nested_monitor.csv):
#   - the transport through the Pacific (north) edge should be about +0.7 to +0.9 Sv into the domain (Bering
#     throughflow), carried out through the Atlantic (south) edge. A sign or size far from that means the velocity
#     rotation or the edge sampling is wrong on this grid;
#   - the domain-mean free surface: plain Flather, no volume constraint, so GLORYS's edge transports that do not
#     close show up as a drift (about 0.5 cm/day per Sv);
#   - the shelf and deep η changes, to compare with our 15 km runs.

include("common.jl")
using CopernicusMarine
using CUDA
using Dates

output_dir = isempty(ARGS) ? joinpath(@__DIR__, "output-1deg") : ARGS[1]
mkpath(output_dir)

grid = panarctic_grid(GPU(); spacing = 100kilometers)

# Only the Atlantic and Pacific edges are open (the default opens every side that touches water). The rest takes the
# branch's defaults: oblique radiation (9 h inflow), tracer reservoirs (20 km), Flather and Chapman, a 2/3° sponge on
# T, S and the baroclinic velocity, the seafloor blended to GLORYS's over 1/3°, and GLORYS interpolated as the initial
# state (not balanced).
dates = (DateTime(1993, 1, 1), DateTime(1993, 1, 6))
ocean = nested_ocean_simulation(grid, GLORYSDaily(); dates, open_sides = (:south, :north),
                                Δt = 20minutes, stop_time = 4days)

child = ocean.model.child
@info "Child grid topology: $(Oceananigans.Grids.topology(child.grid))"

add_nested_monitor!(ocean, child, output_dir)
add_nested_output!(ocean, child, output_dir)

run!(ocean)
