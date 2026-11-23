# Large-eddy simulation of a case in cases.jl: Oceananigans' NonhydrostaticModel, implicit LES with WENO(order = 9)
# (as LESbrary.jl and Oceananigans' Langmuir example), horizontally periodic, a sponge near the bottom, buoyancy
# tracer (or T and S with TEOS-10 for seawater cases). Writes horizontal averages (b or T, S; u, v; and the fluxes
# w b or w T, w S; w u, w v) every 10 minutes to output/<case>_les.jld2. On a GPU, in a REPL from the PanArcticShenanigans clone's environment:
#
#     empty!(ARGS); push!(ARGS, "free_convection")            # or "strong_wind", "strong_wind_stokes"
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/les.jl")
#
# Domain depth from the case (cases.jl: 256 m, 128 m for the Arctic case). A second argument sets the size:
#   (none)  the case's own: 512 × 512 m at 2 m (the calibration LES's extent, 8.4 million cells at 256 m depth), or
#           256 × 256 m at 1 m for arctic_brine_convection (8.4 million cells)
#   small   256 × 256 m at 2 m (2.1 million cells; a few minutes): a quick check
#   fine    256 × 256 m at 1 m (16.8 million cells): when the top metres matter (1 m columns)
# A third argument "cpu" runs on the CPU even when a GPU is there (23 Nov: CPU cost of the LES matrix), e.g.
#     empty!(ARGS); push!(ARGS, "arctic_ice_drag", "small", "cpu")
# Horizontally the domain must stay several times wider than the boundary layer (plumes as wide as it is deep).

using Oceananigans, Oceananigans.Units
using CUDA, JLD2, Printf, Random

include(joinpath(@__DIR__, "cases.jl"))

case_name = isempty(ARGS) ? "free_convection" : ARGS[1]
size_name = length(ARGS) ≥ 2 ? ARGS[2] : "default"
size_name in ("default", "small", "fine") || error("Size must be small or fine (or nothing), not $size_name")
case      = case_parameters(case_name)
arch_name = length(ARGS) ≥ 3 ? ARGS[3] : "auto"     # "cpu" forces the CPU (threads: julia -t N)
arch      = arch_name == "cpu" ? CPU() : CUDA.functional() ? GPU() : CPU()

Δ  = size_name == "fine" ? 1.0 : size_name == "small" ? 2.0 : case.les_spacing   # grid spacing (m)
Lh = size_name == "default" ? case.les_width : 256.0                               # horizontal extent (m)
Lz = case.depth
grid = RectilinearGrid(arch; size = (round(Int, Lh / Δ), round(Int, Lh / Δ), round(Int, Lz / Δ)),
                       x = (0, Lh), y = (0, Lh), z = (-Lz, 0), halo = (5, 5, 5),
                       topology = (Periodic, Periodic, Bounded))

# Surface fluxes (positive upward: Jb > 0 is a buoyancy loss, Js < 0 adds salt, τx < 0 pushes the water towards +x)
momentum_bcs = momentum_boundary_conditions(case)   # constant τx, or quadratic drag on drifting ice (cases.jl)

# Sponge over the bottom tenth, relaxing to the deep (linear) part of the initial profile (as LESbrary)
linear_target(value_at_bottom, gradient) = LinearTarget{:z}(intercept = value_at_bottom + gradient * Lz, gradient = gradient)
mask   = GaussianMask{:z}(center = -Lz, width = Lz / 10)
sponge = Relaxation(rate = 4 / hour, mask = mask)

if case.seawater
    tracer_names = (:T, :S)
    buoyancy_model = SeawaterBuoyancy(equation_of_state = teos10)
    tracer_bcs   = (T = FieldBoundaryConditions(top = FluxBoundaryCondition(case.Jt),
                                                bottom = GradientBoundaryCondition(case.dTdz_bottom)),
                    S = FieldBoundaryConditions(top = FluxBoundaryCondition(case.Js),
                                                bottom = GradientBoundaryCondition(case.dSdz_bottom)))
    tracer_sponges = (T = Relaxation(rate = 4 / hour, mask = mask, target = linear_target(case.T_init(-Lz), case.dTdz_bottom)),
                      S = Relaxation(rate = 4 / hour, mask = mask, target = linear_target(case.S_init(-Lz), case.dSdz_bottom)))
else
    tracer_names = :b
    buoyancy_model = BuoyancyTracer()
    tracer_bcs   = (; b = FieldBoundaryConditions(top = FluxBoundaryCondition(case.Jb),
                                                  bottom = GradientBoundaryCondition(case.N²_bottom)))
    tracer_sponges = (; b = Relaxation(rate = 4 / hour, mask = mask, target = linear_target(case.b_init(-Lz), case.N²_bottom)))
end

# Broadband Stokes drift (cases.jl), its constants as parameters (no globals in GPU kernels; safe to rerun in a REPL)
stokes_drift = nothing
if case.stokes
    stokes = broadband_stokes_parameters(case.τx)
    stokes_drift = UniformStokesDrift(∂z_uˢ = broadband_stokes_shear, parameters = stokes)
    @info @sprintf("Broadband Stokes drift: peak wavelength %.0f m, ∂z uˢ at 1, 10, 30 m: %.2e, %.2e, %.2e s⁻¹",
                   2π / stokes.kᵖ, broadband_stokes_shear(-1.0, 0, stokes), broadband_stokes_shear(-10.0, 0, stokes),
                   broadband_stokes_shear(-30.0, 0, stokes))
end

model = NonhydrostaticModel(grid; coriolis = FPlane(f = case.f),
                            advection = WENO(order = 9),
                            tracers = tracer_names, buoyancy = buoyancy_model,
                            stokes_drift,
                            forcing = (; u = sponge, v = sponge, w = sponge, tracer_sponges...),
                            boundary_conditions = (; momentum_bcs..., tracer_bcs...))

# Initial state at rest, with noise in the top 4 m to start the turbulence (buoyancy noise 2.56e-5 m s⁻², as the
# three-layer cases had: 1e-2 × N²_thermocline × 256 m; for seawater as salinity noise of 3e-3 g/kg ≈ 2.2e-5 m s⁻²)
Random.seed!(1234)
noise(z) = randn() * exp(z / 4)
u★ = max(sqrt(abs(case.τx)), sqrt(case.Cd) * hypot(case.uᵢ, case.vᵢ), 1e-3)
if case.seawater
    set!(model, T = (x, y, z) -> case.T_init(z), S = (x, y, z) -> case.S_init(z) + 3e-3 * noise(z),
                u = (x, y, z) -> 1e-1 * u★ * noise(z), w = (x, y, z) -> 1e-1 * u★ * noise(z))
else
    set!(model, b = (x, y, z) -> case.b_init(z) + 2.56e-5 * noise(z),
                u = (x, y, z) -> 1e-1 * u★ * noise(z), w = (x, y, z) -> 1e-1 * u★ * noise(z))
end

simulation = Simulation(model; Δt = 1.0, stop_time = case.stop_time)
conjure_time_step_wizard!(simulation; cfl = 0.7, max_Δt = 30.0)

wall = Ref(time_ns())
function progress(sim)
    u, v, w = sim.model.velocities
    @info @sprintf("%s: t = %s, Δt = %s, max|w| = %.2e m/s, wall %s", case_name, prettytime(sim), prettytime(sim.Δt),
                   maximum(abs, w), prettytime(1e-9 * (time_ns() - wall[])))
    return nothing
end
add_callback!(simulation, progress, IterationInterval(500))

u, v, w = model.velocities
momentum_averages = (; u = Average(u, dims = (1, 2)), v = Average(v, dims = (1, 2)),
                       wu = Average(w * u, dims = (1, 2)), wv = Average(w * v, dims = (1, 2)))
if case.seawater
    T, S = model.tracers.T, model.tracers.S
    tracer_averages = (; T = Average(T, dims = (1, 2)), S = Average(S, dims = (1, 2)),
                         wT = Average(w * T, dims = (1, 2)), wS = Average(w * S, dims = (1, 2)))
else
    b = model.tracers.b
    tracer_averages = (; b = Average(b, dims = (1, 2)), wb = Average(w * b, dims = (1, 2)))
end
averages = (; tracer_averages..., momentum_averages...)
file = joinpath(output_folder(), case_name * (size_name == "default" ? "" : "_" * size_name) * "_les.jld2")
simulation.output_writers[:averages] = JLD2Writer(model, averages; filename = file,
                                                  schedule = TimeInterval(10minutes), overwrite_files = true)

@info "LES $case_name on $(arch isa GPU ? "GPU" : "CPU"): $(size(grid)) cells of $(Δ) m; output $file"
run!(simulation)
@info "Done: $file"
