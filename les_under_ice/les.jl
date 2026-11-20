# Large-eddy simulation of a case in cases.jl: Oceananigans' NonhydrostaticModel, implicit LES with WENO(order = 9)
# (as LESbrary.jl and Oceananigans' Langmuir example), horizontally periodic, a sponge near the bottom, buoyancy
# tracer. Writes horizontal averages (b, u, v, and the fluxes w b, w u, w v) every 10 minutes to
# output/<case>_les.jld2. On a GPU, in a REPL from the PanArcticShenanigans clone's environment:
#
#     empty!(ARGS); push!(ARGS, "free_convection")            # or "strong_wind", "strong_wind_stokes"
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/les.jl")
#
# Domain 256 m deep (cases.jl). A second argument sets the size:
#   (none)  512 × 512 m at 2 m: the calibration LES's extent (8.4 million cells; ~15–30 min on an A100 for 24 hours)
#   small   256 × 256 m at 2 m (2.1 million cells; a few minutes): a quick check
#   fine    256 × 256 m at 1 m (16.8 million cells): when the top metres matter (1 m columns)
# Horizontally the domain must stay several times wider than the boundary layer (plumes as wide as it is deep).

using Oceananigans, Oceananigans.Units
using CUDA, JLD2, Printf, Random

include(joinpath(@__DIR__, "cases.jl"))

case_name = isempty(ARGS) ? "free_convection" : ARGS[1]
size_name = length(ARGS) ≥ 2 ? ARGS[2] : "default"
size_name in ("default", "small", "fine") || error("Size must be small or fine (or nothing), not $size_name")
case      = case_parameters(case_name)
arch      = CUDA.functional() ? GPU() : CPU()

Δ  = size_name == "fine" ? 1.0 : 2.0           # grid spacing (m)
Lh = size_name == "default" ? 512.0 : 256.0    # horizontal extent (m)
Lz = domain_depth
grid = RectilinearGrid(arch; size = (round(Int, Lh / Δ), round(Int, Lh / Δ), round(Int, Lz / Δ)),
                       x = (0, Lh), y = (0, Lh), z = (-Lz, 0), halo = (5, 5, 5),
                       topology = (Periodic, Periodic, Bounded))

# Surface fluxes (positive upward: Jb > 0 is a buoyancy loss, τx < 0 pushes the water towards +x)
b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(case.Jb), bottom = GradientBoundaryCondition(three_layer.N²_deep))
u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(case.τx))

# Sponge over the bottom tenth, relaxing to the deep stratification (as LESbrary)
z_deep   = -three_layer.surface_layer_depth - three_layer.thermocline_width
b_target = LinearTarget{:z}(intercept = three_layer_buoyancy(z_deep) - three_layer.N²_deep * z_deep,
                            gradient  = three_layer.N²_deep)
mask     = GaussianMask{:z}(center = -Lz, width = Lz / 10)
sponge   = Relaxation(rate = 4 / hour, mask = mask)
b_sponge = Relaxation(rate = 4 / hour, mask = mask, target = b_target)

# Stokes drift shear, with its constants as parameters (no globals in GPU kernels; safe to rerun in a REPL)
@inline monochromatic_stokes_shear(z, t, p) = 2 * p.k * p.uˢ₀ * exp(2 * p.k * z)
stokes_drift = nothing
if case.stokes
    stokes = (uˢ₀ = stokes_surface_velocity(case.τx), k = stokes_wavenumber)
    stokes_drift = UniformStokesDrift(∂z_uˢ = monochromatic_stokes_shear, parameters = stokes)
    @info @sprintf("Stokes drift: uˢ₀ = %.3f m/s, decay scale %.1f m (La_t = 0.3)", stokes.uˢ₀, 1 / 2stokes.k)
end

model = NonhydrostaticModel(grid; coriolis = FPlane(f = case.f),
                            advection = WENO(order = 9),
                            tracers = :b, buoyancy = BuoyancyTracer(),
                            stokes_drift,
                            forcing = (u = sponge, v = sponge, w = sponge, b = b_sponge),
                            boundary_conditions = (u = u_bcs, b = b_bcs))

# Initial state: three-layer stratification at rest, with noise in the top 4 m to start the turbulence
Random.seed!(1234)
noise(z) = randn() * exp(z / 4)
u★ = max(sqrt(abs(case.τx)), 1e-3)
set!(model, b = (x, y, z) -> three_layer_buoyancy(z) + 1e-2 * three_layer.N²_thermocline * Lz * noise(z),
            u = (x, y, z) -> 1e-1 * u★ * noise(z), w = (x, y, z) -> 1e-1 * u★ * noise(z))

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
b = model.tracers.b
averages = (; b = Average(b, dims = (1, 2)), u = Average(u, dims = (1, 2)), v = Average(v, dims = (1, 2)),
              wb = Average(w * b, dims = (1, 2)), wu = Average(w * u, dims = (1, 2)), wv = Average(w * v, dims = (1, 2)))
file = joinpath(output_folder(), case_name * (size_name == "default" ? "" : "_" * size_name) * "_les.jld2")
simulation.output_writers[:averages] = JLD2Writer(model, averages; filename = file,
                                                  schedule = TimeInterval(10minutes), overwrite_files = true)

@info "LES $case_name on $(arch isa GPU ? "GPU" : "CPU"): $(size(grid)) cells of $(Δ) m; output $file"
run!(simulation)
@info "Done: $file"
