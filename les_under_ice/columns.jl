# Single-column simulations of a case in cases.jl: Oceananigans' HydrostaticFreeSurfaceModel on a (Flat, Flat,
# Bounded) grid with CATKE, the same forcing and initial state as les.jl, at 1, 4 and 8 m and on the pan-Arctic model's
# vertical grids at 50 and 100 levels (top cells ~8 and ~4 m). CATKE variants:
#   stock            CATKE as calibrated (surface TKE flux Cᵂu★ u★³ + CᵂwΔ (wΔ)³)
#   no_surface_shear Cᵂu★ = 0: no shear-driven (wave-like) surface TKE flux; the convective one is kept
#   rigid_wall       Cᵂu★ = 0 and Cˢ = 0.967 (von Kármán constant 0.4, as at the bottom): CATKE for a wall without waves
#   k_epsilon        Oceananigans' two-equation k-ε closure (TKEDissipationVerticalDiffusivity, default constants;
#                    whether they come from a calibration against LESbrary is not known here): an independent model
# Writes b (or T and S), u, v, e and κᶜ every 10 minutes to output/<case>_column_<variant>_<dz1|…|nz50|nz100>.jld2. Fast on a CPU:
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
# Uniform 1, 4, 8 m, and the pan-Arctic model's own grids at 50 and 100 levels (cases.jl, column_faces)
spacings = column_spacings

# CATKE adds its TKE tracer e itself. In a function, so a REPL rerun has no soft-scope warnings.
function run_column(case_name, case, variant, closure, spacing)
    grid = RectilinearGrid(CPU(); size = length(column_faces(spacing, case.depth)) - 1,
                           z = column_faces(spacing, case.depth), topology = (Flat, Flat, Bounded))
    momentum_bcs = momentum_boundary_conditions(case; column = true)
    if case.seawater
        tracers  = (:T, :S)
        buoyancy = SeawaterBuoyancy(equation_of_state = teos10)
        tracer_bcs = (T = FieldBoundaryConditions(top = FluxBoundaryCondition(case.Jt),
                                                  bottom = GradientBoundaryCondition(case.dTdz_bottom)),
                      S = FieldBoundaryConditions(top = FluxBoundaryCondition(case.Js),
                                                  bottom = GradientBoundaryCondition(case.dSdz_bottom)))
    else
        tracers  = :b
        buoyancy = BuoyancyTracer()
        tracer_bcs = (; b = FieldBoundaryConditions(top = FluxBoundaryCondition(case.Jb),
                                                    bottom = GradientBoundaryCondition(case.N²_bottom)))
    end
    model = HydrostaticFreeSurfaceModel(grid; closure, coriolis = FPlane(f = case.f), tracers, buoyancy,
                                        boundary_conditions = (; momentum_bcs..., tracer_bcs...))
    if case.seawater
        set!(model, T = case.T_init, S = case.S_init, e = 1e-6)
        tracer_outputs = (; T = model.tracers.T, S = model.tracers.S)
    else
        set!(model, b = case.b_init, e = 1e-6)
        tracer_outputs = (; b = model.tracers.b)
    end
    simulation = Simulation(model; Δt = 1minute, stop_time = case.stop_time)
    file = joinpath(output_folder(), @sprintf("%s_column_%s_%s.jld2", case_name, variant, spacing_label(spacing)))
    outputs = (; tracer_outputs..., u = model.velocities.u, v = model.velocities.v, e = model.tracers.e,
                 κc = model.closure_fields.κc)
    simulation.output_writers[:profiles] = JLD2Writer(model, outputs; filename = file,
                                                      schedule = TimeInterval(10minutes), overwrite_files = true)
    run!(simulation)
    @info "Column $case_name, $variant, $(spacing_label(spacing)): $file"
    return file
end

for (variant, closure) in sort(collect(variants); by = first), spacing in spacings
    run_column(case_name, case, variant, closure, spacing)
end
