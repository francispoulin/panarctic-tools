# Why do the fixes not take effect? Builds the f(72°N) partial-cell ridge model once and checks, without time
# stepping, whether each override is reached and what it changes.
#
#   julia --project=<PanArcticShenanigans>/nested probe_fixes.jl [profile.csv] [partial|whole]

using Oceananigans, Printf
using Oceananigans.Architectures: CPU, convert_to_device   # the conversion launch! applies to kernel arguments
using Oceananigans.Operators: Operators
using Oceananigans.Grids: inactive_node, Center

RIDGE_NO_RUN = true
include(joinpath(@__DIR__, "ridge_fixes.jl"))   # definitions only: ridge_model, the fixes

bottom_kind = length(ARGS) ≥ 2 ? Symbol(ARGS[2]) : :partial
model = ridge_model(bottom_kind, 20, 72)
time_step!(model, 10minutes)   # so that the grid stretching σ (and ∂x z) reflects η; at t = 0 it may still be 1
grid  = model.grid
pHY′  = model.pressure.pHY′

kgrid, kp = convert_to_device(CPU(), (grid, pHY′))
println("bottom: ", bottom_kind)

println("\n--- types ---")
println("model.grid isa BottomMutableGrid:         ", grid  isa BottomMutableGrid)
println("adapted grid isa BottomMutableGrid:       ", kgrid isa BottomMutableGrid)
println("adapted pHY′ isa OffsetArray3:            ", kp    isa OffsetArray3)
println("typeof(adapted pHY′):                     ", typeof(kp))
println("z topology:                               ", Oceananigans.Grids.topology(grid, 3))
println("typeof(grid.underlying_grid.z):           ", typeof(grid.underlying_grid.z).name.name)
println("\nmethod ∂xᶠᶜᶜ picks for the kernel arguments:")
println("  ", which(Operators.∂xᶠᶜᶜ, Tuple{Int, Int, Int, typeof(kgrid), typeof(kp)}))
println("method update_hydrostatic_pressure! picks:")
println("  ", which(update_hydrostatic_pressure!, Tuple{typeof(pHY′), typeof(model.architecture), typeof(grid),
                                                     typeof(model.buoyancy), typeof(model.tracers)}))

# Force at the cliff faces with each switch, at t = 0 (η tilted across the channel).
function cliff_force(i)
    out = Float64[]
    for k in 1:size(grid, 3), j in 1:size(grid, 2)
        push!(out, Operators.∂xᶠᶜᶜ(i, j, k, kgrid, kp))
    end
    return out
end

function pressure_with(; F1 = false, F2 = false)
    set_fixes!(; F1, F2)
    update_hydrostatic_pressure!(pHY′, model.architecture, grid, model.buoyancy, model.tracers)
    Oceananigans.BoundaryConditions.fill_halo_regions!(pHY′)
    return copy(interior(pHY′)), vcat(cliff_force(13), cliff_force(21))
end

ζ = [Operators.∂x_zᶠᶜᶜ(i, j, k, kgrid) for i in (13, 21), j in 1:size(grid, 2), k in 1:size(grid, 3)]
@printf("\nmax |∂x z| at the cliff faces:              %.3e\n", maximum(abs, ζ))

p0, f0 = pressure_with()
p2, f2 = pressure_with(F2 = true)
p1, f1 = pressure_with(F1 = true)
calls_F1 = F1_CALLS[]
p12, f12 = pressure_with(F1 = true, F2 = true)
set_fixes!()

println("\n--- effect at t = 0 ---")
@printf("max |Δp′| from F2 (all cells):              %.3e m²/s²\n", maximum(abs, p2 .- p0))
@printf("max |∂x p′| at the cliff faces, stock:      %.3e m/s²\n", maximum(abs, f0))
@printf("max |Δ force| at cliff faces, F1:           %.3e m/s²\n", maximum(abs, f1 .- f0))
@printf("max |Δ force| at cliff faces, F2:           %.3e m/s²\n", maximum(abs, f2 .- f0))
@printf("max |Δ force| at cliff faces, F1+F2:        %.3e m/s²\n", maximum(abs, f12 .- f0))
println("F1 override calls during the F1 probe: ", calls_F1, "  (", 2 * size(grid, 2) * size(grid, 3), " expected)")
