# run_column: one single column of a case (cases.jl) with a given closure on a given grid, writing b (or T, S), u, v, e
# and the diffusivities κc, κu every 10 minutes. Used by columns.jl and variants.jl.

using Oceananigans, Oceananigans.Units
using JLD2, Printf

# The model of one column: grid (uniform spacing or a pan-Arctic grid), forcing and initial state of the case
function column_model(case, closure, spacing)
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
    case.seawater ? set!(model, T = case.T_init, S = case.S_init, e = 1e-6) : set!(model, b = case.b_init, e = 1e-6)
    return model
end

# The buoyancy profile at the end of a run (b, or from T and S for seawater), without writing output: for screen.jl
function final_buoyancy(case, closure, spacing; Δt = 1minute)
    model = column_model(case, closure, spacing)
    simulation = Simulation(model; Δt, stop_time = case.stop_time, verbose = false)
    run!(simulation)
    if case.seawater
        T, S = vec(Array(interior(model.tracers.T))), vec(Array(interior(model.tracers.S)))
        return seawater_buoyancy.(T, S)
    else
        return vec(Array(interior(model.tracers.b)))
    end
end

# CATKE adds its TKE tracer e itself. In a function, so a REPL rerun has no soft-scope warnings.
function run_column(case_name, case, variant, closure, spacing; Δt = 1minute,
                    file = joinpath(output_folder(), @sprintf("%s_column_%s_%s.jld2", case_name, variant, spacing_label(spacing))))
    model = column_model(case, closure, spacing)
    tracer_outputs = case.seawater ? (; T = model.tracers.T, S = model.tracers.S) : (; b = model.tracers.b)
    simulation = Simulation(model; Δt, stop_time = case.stop_time)
    # Diffusivities on faces: κc for tracers, κu for momentum (both CATKE and k-ε have them), for diagnose_fluxes.jl
    diffusivities = (; κc = model.closure_fields.κc,
                       (hasproperty(model.closure_fields, :κu) ? (; κu = model.closure_fields.κu) : (;))...)
    outputs = (; tracer_outputs..., u = model.velocities.u, v = model.velocities.v, e = model.tracers.e, diffusivities...)
    simulation.output_writers[:profiles] = JLD2Writer(model, outputs; filename = file,
                                                      schedule = TimeInterval(10minutes), overwrite_files = true)
    run!(simulation)
    @info "Column $case_name, $variant, $(spacing_label(spacing)): $file"
    return file
end

