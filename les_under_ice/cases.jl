# Cases for the LES-vs-single-column comparison (claude/catke-les-plan.md). Included by les.jl, columns.jl, compare.jl.
#
# Step 0 (this file, for now): CATKE's own calibration regime, to check the pipeline before adding ice. Two cases of
# the 24-hour suite of LESbrary.jl (src/IdealizedExperiments/IdealizedExperiments.jl), used by Wagner et al. (arXiv
# 2306.13204) to calibrate CATKE, with LESbrary's "three-layer" initial stratification and a buoyancy tracer:
#   free_convection     surface buoyancy loss 2.4e-7 m² s⁻³, no wind
#   strong_wind         wind stress 6.8e-4 m² s⁻² (u★ ≈ 0.026 m/s), no buoyancy flux, no waves
#   strong_wind_stokes  the same with surface-wave Stokes drift (Langmuir turbulence), as in the calibration
# All with f = 1e-4 s⁻¹, 24 hours. LESbrary's calibration LES had a broadband Stokes drift (ConstantFluxStokesDrift);
# here a monochromatic wave with turbulent Langmuir number La_t = (u★ / uˢ₀)^½ = 0.3, as in Oceananigans' Langmuir
# example: close to, not identical with, the calibration. strong_wind against strong_wind_stokes shows how much of
# CATKE's surface TKE flux (Cᵂu★ u★³) stands for waves: the question under sea ice.
#
# Stage 1 (sea ice, GLORYS profile, TEOS-10, ice–ocean drag) will be added here after step 0 checks out.

using Oceananigans.Units

three_layer = (; surface_layer_depth = 48.0, thermocline_width = 24.0,
                 N²_surface = 2e-6, N²_thermocline = 1e-5, N²_deep = 2e-6)

function three_layer_buoyancy(z; p = three_layer)
    h, w = p.surface_layer_depth, p.thermocline_width
    z ≥ -h     && return p.N²_surface * z
    z ≥ -h - w && return -p.N²_surface * h + p.N²_thermocline * (z + h)
    return -p.N²_surface * h - p.N²_thermocline * w + p.N²_deep * (z + h + w)
end

cases = Dict(
    "free_convection"    => (; τx = 0.0,     Jb = 2.4e-7, f = 1e-4, stop_time = 24hours, stokes = false),
    "strong_wind"        => (; τx = -6.8e-4, Jb = 0.0,    f = 1e-4, stop_time = 24hours, stokes = false),
    "strong_wind_stokes" => (; τx = -6.8e-4, Jb = 0.0,    f = 1e-4, stop_time = 24hours, stokes = true))

# Monochromatic Stokes drift with La_t = 0.3: uˢ(z) = uˢ₀ exp(2kz), wavelength 60 m (decay scale 1/2k ≈ 4.8 m)
stokes_wavenumber = 2π / 60
stokes_surface_velocity(τx; La_t = 0.3) = sqrt(abs(τx)) / La_t^2

domain_depth = 128.0   # m: enough for 24 hours of these forcings (the three-layer thermocline is at 48–72 m)

function case_parameters(name)
    haskey(cases, name) || error("Unknown case \"$name\"; cases: $(join(sort(collect(keys(cases))), ", "))")
    return cases[name]
end

output_folder() = mkpath(joinpath(@__DIR__, "output"))
