# Single-column simulations of a case in cases.jl: Oceananigans' HydrostaticFreeSurfaceModel on a (Flat, Flat,
# Bounded) grid with CATKE, the same forcing and initial state as les.jl, at 1, 4 and 8 m (8 m is the pan-Arctic
# model's top cell at 50 levels). CATKE variants:
#   stock            CATKE as calibrated (surface TKE flux Cᵂu★ u★³ + CᵂwΔ (wΔ)³)
#   no_surface_shear Cᵂu★ = 0: no shear-driven (wave-like) surface TKE flux; the convective one is kept
#   rigid_wall       Cᵂu★ = 0 and Cˢ = 0.967 (von Kármán constant 0.4, as at the bottom): CATKE for a wall without waves
#   k_epsilon        Oceananigans' two-equation k-ε closure (TKEDissipationVerticalDiffusivity, default constants;
#                    whether they come from a calibration against LESbrary is not known here): an independent model
# Writes b, u, v, e and κᶜ every 10 minutes to output/<case>_column_<variant>_dz<Δz>.jld2. Fast on a CPU:
#
#     empty!(ARGS); push!(ARGS, "free_convection")
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/columns.jl")

using Oceananigans, Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity, TKEDissipationVerticalDiffusivity
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEEquation, CATKEMixingLength
using JLD2, Printf

include(joinpath(@__DIR__, "cases.jl"))

case_name = isempty(ARGS) ? "free_convection" : ARGS[1]
case      = case_parameters(case_name)

# rigid_wall: CATKE's surface mixing-length coefficient Cˢ = 1.131 corresponds to a von Kármán constant κ = 0.47,
# larger than a rigid wall's 0.4 "because the large eddy simulations used to calibrate CATKE include surface wave
# effects" (Oceananigans, catke_mixing_length.jl); its bottom coefficient Cᵇ = 0.967 gives κ = 0.4. Under ice (a rigid
# wall, no waves): Cˢ = Cᵇ = 0.967 and no shear-driven surface TKE flux.
variants = Dict("stock"            => CATKEVerticalDiffusivity(),
                "no_surface_shear" => CATKEVerticalDiffusivity(turbulent_kinetic_energy_equation = CATKEEquation(Cᵂu★ = 0.0)),
                "rigid_wall"       => CATKEVerticalDiffusivity(mixing_length = CATKEMixingLength(Cˢ = 0.967),
                                                               turbulent_kinetic_energy_equation = CATKEEquation(Cᵂu★ = 0.0)),
                "k_epsilon"        => TKEDissipationVerticalDiffusivity())
spacings = (1, 4, 8)

# CATKE adds its TKE tracer e itself (tracers = :b only). In a function, so a REPL rerun has no soft-scope warnings.
function run_column(case_name, case, variant, closure, Δz)
    grid = RectilinearGrid(CPU(); size = round(Int, domain_depth / Δz), z = (-domain_depth, 0),
                           topology = (Flat, Flat, Bounded))
    b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(case.Jb),
                                    bottom = GradientBoundaryCondition(three_layer.N²_deep))
    u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(case.τx))
    model = HydrostaticFreeSurfaceModel(grid; closure, coriolis = FPlane(f = case.f),
                                        tracers = :b, buoyancy = BuoyancyTracer(),
                                        boundary_conditions = (u = u_bcs, b = b_bcs))
    set!(model, b = three_layer_buoyancy, e = 1e-6)
    simulation = Simulation(model; Δt = 1minute, stop_time = case.stop_time)
    file = joinpath(output_folder(), @sprintf("%s_column_%s_dz%d.jld2", case_name, variant, Δz))
    outputs = (; b = model.tracers.b, u = model.velocities.u, v = model.velocities.v, e = model.tracers.e,
                 κc = model.closure_fields.κc)
    simulation.output_writers[:profiles] = JLD2Writer(model, outputs; filename = file,
                                                      schedule = TimeInterval(10minutes), overwrite_files = true)
    run!(simulation)
    @info "Column $case_name, $variant, Δz = $Δz m: $file"
    return file
end

for (variant, closure) in sort(collect(variants); by = first), Δz in spacings
    run_column(case_name, case, variant, closure, Δz)
end
