# Cases for the LES-vs-single-column comparison (claude/catke-les-plan.md). Included by les.jl, columns.jl, compare.jl.
#
# Step 0 (this file, for now): CATKE's own calibration regime, to check the pipeline before adding ice. Two cases of
# the 24-hour suite of LESbrary.jl (src/IdealizedExperiments/IdealizedExperiments.jl), used by Wagner et al. (arXiv
# 2306.13204) to calibrate CATKE, with LESbrary's "three-layer" initial stratification and a buoyancy tracer:
#   free_convection     surface buoyancy loss 2.4e-7 m² s⁻³, no wind
#   strong_wind         wind stress 6.8e-4 m² s⁻² (u★ ≈ 0.026 m/s), no buoyancy flux, no waves
#   strong_wind_stokes  the same with the calibration's broadband surface-wave Stokes drift (Langmuir turbulence);
#                       this is LESbrary's strong_wind (step 1: compared with their LES in compare_lesbrary.jl)
# All with f = 1e-4 s⁻¹, 24 hours. strong_wind against strong_wind_stokes shows how much of CATKE's mixing stands for
# waves: CATKE was calibrated only on cases with Stokes drift (every wind case in LESbrary has it).
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

# Broadband Stokes drift of the calibration LES (LESbrary.jl's ConstantFluxStokesDrift: the Lenain and Pizzo 2020
# spectrum for wind seas in equilibrium with a constant wind stress, peak wavenumber kᵖ = 1e-6 g / |τ|). Step 1 (21 Nov)
# showed that the monochromatic wave first used here (decay scale 4.8 m) gave a boundary layer 10 m shallower than
# their LES: their Stokes drift reaches tens of metres (0.074 m/s at 10 m, 0.038 at 20 m, for the strong_wind stress).
# The shear ∂z uˢ is written in closed form (derived from their uˢ(z); checked against a numerical derivative), with
# an erf approximation (Abramowitz and Stegun 7.1.26, error < 2e-7), so it runs on the GPU with no table.
function broadband_stokes_parameters(τx; g = 9.81, water_density = 1024.0, air_density = 1.225,
                                     Cᵝ = 0.105, Cʳ = 9.7e-3, Cⁱ = 0.072, Cᴮ = 7e-3)
    a★ = sqrt(water_density * abs(τx) / air_density)   # air friction velocity
    return (A = a★ * Cᵝ, B = 2Cᴮ * sqrt(g),
            kᵖ = 1e-6 * g / abs(τx), kⁿ = Cʳ * g / a★^2, kⁱ = Cⁱ * g / a★^2, z_min = 0.05)
end

@inline function erf_approximation(x)
    t = 1 / (1 + 0.3275911x)
    return 1 - t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) * exp(-x^2)
end

# ∂z uˢ for uˢ = A [E₁ terms between kᵖ and kⁿ] + B [saturation range between kⁿ and kⁱ] (Lenain and Pizzo 2020)
@inline function broadband_stokes_shear(z, t, p)
    z = min(z, -p.z_min)   # finite at the surface face (where it multiplies w = 0)
    s = -z
    equilibrium = p.A * (exp(2 * p.kⁿ * z) - exp(2 * p.kᵖ * z)) / z
    saturation  = p.B * sqrt(2π / s) * (erf_approximation(sqrt(2 * p.kⁱ * s)) - erf_approximation(sqrt(2 * p.kⁿ * s)))
    return equilibrium + saturation
end

# 256 m deep, as the calibration LES: free convection deepens to ~100 m in 24 hours through the weak deep stratification,
# and the bottom sponge (lowest tenth) must stay below the boundary layer (step 0 at 128 m: it did not)
domain_depth = 256.0

function case_parameters(name)
    haskey(cases, name) || error("Unknown case \"$name\"; cases: $(join(sort(collect(keys(cases))), ", "))")
    return cases[name]
end

output_folder() = mkpath(joinpath(@__DIR__, "output"))
