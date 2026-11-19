# The thin-bottom-cell cases with plain Oceananigans, no overrides: for testing an Oceananigans branch as it is, such
# as #6538 (vector-invariant divergence with the area difference interpolated symmetrically) against its base.
#
#   julia --project=<environment with the branch> branch_checks.jl [days]
#
# Four cases, each whole vs partial cells, WENOVectorInvariant() ("stock", i.e. the branch's own) and the flux form
# (reference): the 2-D channel of #6532 without rotation and with f at 72°N, the same channel along y at 72°N (the
# v-face code), and the 3-D seamount at 72°N. Same numbers as the #6533 PR table for "stock" on main: 2.39 / 0.40 m/s,
# 1.15 / 0.37, and 2.11 / 0.85 for the seamount. CPU; the seamount takes some minutes per run.

const PLAIN_OCEANANIGANS = true
include(joinpath(@__DIR__, "thin_cell_cases.jl"))

days = length(ARGS) ≥ 1 ? parse(Float64, ARGS[1]) : 10.0
cases = [ThinCellCase(geometry = :x),
         ThinCellCase(geometry = :x, latitude = 72),
         ThinCellCase(geometry = :y, latitude = 72),
         ThinCellCase(geometry = :xy, latitude = 72)]

for c in cases
    scheme_table(c, ["stock", "flux"]; days)
end
