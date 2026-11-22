# Cases for the LES-vs-single-column comparison (claude/catke-les-plan.md). Included by les.jl, columns.jl, compare.jl.
#
# Steps 0–2: CATKE's own calibration regime, before adding ice. The 24-hour suite of LESbrary.jl
# (src/IdealizedExperiments/IdealizedExperiments.jl), used by Wagner et al. (arXiv 2306.13204) to calibrate CATKE, with
# LESbrary's "three-layer" initial stratification and a buoyancy tracer (rotating cases only; free_convection has no
# wind, so no waves), e.g.:
#   free_convection     surface buoyancy loss 2.4e-7 m² s⁻³, no wind
#   strong_wind         wind stress 6.8e-4 m² s⁻² (u★ ≈ 0.026 m/s), no buoyancy flux, no waves
#   strong_wind_stokes  the same with the calibration's broadband surface-wave Stokes drift (Langmuir turbulence);
#                       this is LESbrary's strong_wind (step 1: compared with their LES in compare_lesbrary.jl)
# All with f = 1e-4 s⁻¹, 24 hours. strong_wind against strong_wind_stokes shows how much of CATKE's mixing stands for
# waves: CATKE was calibrated only on cases with Stokes drift (every wind case in LESbrary has it).
#
# Stage 1 (sea ice, GLORYS profile, TEOS-10, ice–ocean drag) will be added here after step 0 checks out.

using Oceananigans, Oceananigans.Units
using Oceananigans.OutputReaders: FieldTimeSeries
using JLD2

three_layer = (; surface_layer_depth = 48.0, thermocline_width = 24.0,
                 N²_surface = 2e-6, N²_thermocline = 1e-5, N²_deep = 2e-6)

function three_layer_buoyancy(z; p = three_layer)
    h, w = p.surface_layer_depth, p.thermocline_width
    z ≥ -h     && return p.N²_surface * z
    z ≥ -h - w && return -p.N²_surface * h + p.N²_thermocline * (z + h)
    return -p.N²_surface * h - p.N²_thermocline * w + p.N²_deep * (z + h + w)
end

# Every case has the same fields (helpers below), so the scripts never need to know which kind it is:
#   τx, Jb, Js, Jt      surface fluxes, positive upward: momentum (τx < 0: wind towards +x), buoyancy (Jb > 0: loss),
#                       salt (Js < 0: brine rejection adds salt), heat as temperature flux (Jt > 0: cooling)
#   f, stop_time, stokes
#   seawater            false: buoyancy tracer b with initial b_init(z); true: T and S with TEOS-10 (T_init, S_init)
#   N²_bottom, dTdz_bottom, dSdz_bottom   bottom gradients (boundary condition and the LES sponge's linear target)
#   depth               domain depth (m)
#   les_spacing, les_width   the LES grid spacing and horizontal extent (m) when les.jl is given no size
function buoyancy_case(; τx = 0.0, Jb = 0.0, f = 1e-4, stop_time = 24hours, stokes = false,
                         b_init = three_layer_buoyancy, N²_bottom = three_layer.N²_deep, depth = 256.0,
                         les_spacing = 2.0, les_width = 512.0)
    return (; τx, Jb, Js = 0.0, Jt = 0.0, f, stop_time, stokes, seawater = false,
              b_init, T_init = z -> 0.0, S_init = z -> 0.0, N²_bottom, dTdz_bottom = 0.0, dSdz_bottom = 0.0,
              depth, les_spacing, les_width)
end

function seawater_case(; τx = 0.0, Js = 0.0, Jt = 0.0, f = 1.4e-4, stop_time = 4days, stokes = false,
                         T_init, S_init, dTdz_bottom, dSdz_bottom, depth = 128.0,
                         les_spacing = 1.0, les_width = 256.0)
    return (; τx, Jb = 0.0, Js, Jt, f, stop_time, stokes, seawater = true,
              b_init = z -> 0.0, T_init, S_init, N²_bottom = 0.0, dTdz_bottom, dSdz_bottom,
              depth, les_spacing, les_width)
end

cases = Dict{String, NamedTuple}()

# LESbrary's 24-hour suite (IdealizedExperiments.jl), each wind case without waves (the situation under ice) and with
# the calibration's broadband Stokes drift (suffix _stokes; this is LESbrary's own case, compared in
# compare_lesbrary.jl). Three-layer stratification, f = 1e-4 s⁻¹, 24 hours.
twenty_four_hour_suite = Dict(
    "free_convection"          => (τx = 0.0,     Jb = 2.4e-7),
    "weak_wind_strong_cooling" => (τx = -3.0e-4, Jb = 2.0e-7),
    "med_wind_med_cooling"     => (τx = -4.5e-4, Jb = 1.6e-7),
    "strong_wind_weak_cooling" => (τx = -5.9e-4, Jb = 1.0e-7),
    "strong_wind"              => (τx = -6.8e-4, Jb = 0.0))

for (name, forcing) in twenty_four_hour_suite
    cases[name] = buoyancy_case(; forcing...)
    forcing.τx == 0 || (cases[name * "_stokes"] = buoyancy_case(; forcing..., stokes = true))
end

# Rung 1 of the validation ladder (claude/catke-les-plan.md): free convection into uniform stratification, to test the
# entrainment law h² N² = 2(1 + 2A) Jb t (A ≈ 0.2: the entrainment flux at the base is −A Jb) over a range of forcings,
# not only one case. Chosen so that h(24 h) ≈ 55–110 m (well above the sponge): (N², Jb)
convection_suite = [(1e-5, 2.4e-7), (1e-5, 4.8e-7), (1e-5, 1.2e-7), (4e-5, 4.8e-7), (2.5e-6, 1.2e-7)]
convection_case_name(N², Jb) = "convection_N2_$(N²)_Jb_$(Jb)"
for (N², Jb) in convection_suite
    cases[convection_case_name(N², Jb)] = buoyancy_case(; Jb, b_init = z -> N² * z, N²_bottom = N²)
end

# Rung 1, Arctic version: brine-driven convection under ice at freeze-up, no wind (ice at rest), TEOS-10, f at 75° N.
# A weakly stratified winter mixed layer (N² ≈ 1e-5) over a strong halocline (N² ≈ 2e-4), at the freezing point
# (T = −1.6 °C, Conservative Temperature; Absolute Salinity). Salt flux of ice growth at 5 cm/day:
# Js = −growth × (S_ocean − S_ice) × ρ_ice / ρ_water ≈ −5.8e-7 m/s × 25 g/kg × 0.9 ≈ −1.3e-5 g/kg m/s
# (buoyancy flux ≈ 1e-7 m² s⁻³). No heat flux: the latent heat of growth balances the heat lost to the atmosphere.
# The mixed layer convects to the halocline in about a day; the question is how much the halocline is then eroded
# (the regime where CATKE over-mixes at freeze-up in the pan-Arctic runs). 4 days, 128 m deep, LES at 1 m.
arctic_profile = (; mixed_layer_depth = 40.0, S_surface = 30.0, dSdz_mixed_layer = 1.3e-3, dSdz_halocline = 0.025,
                    T = -1.6)
arctic_salinity(z; p = arctic_profile) =
    -z ≤ p.mixed_layer_depth ? p.S_surface - p.dSdz_mixed_layer * z :
    p.S_surface + p.dSdz_mixed_layer * p.mixed_layer_depth + p.dSdz_halocline * (-z - p.mixed_layer_depth)
cases["arctic_brine_convection"] = seawater_case(; Js = -1.3e-5, T_init = z -> arctic_profile.T,
                                                   S_init = arctic_salinity, dTdz_bottom = 0.0,
                                                   dSdz_bottom = -arctic_profile.dSdz_halocline)

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
# and the bottom sponge (lowest tenth) must stay below the boundary layer (step 0 at 128 m: it did not). Each case has
# its own depth (case.depth); this global is kept for scan.jl.
domain_depth = 256.0

function case_parameters(name)
    haskey(cases, name) || error("Unknown case \"$name\"; cases: $(join(sort(collect(keys(cases))), ", "))")
    return cases[name]
end

output_folder() = mkpath(joinpath(@__DIR__, "output"))

# TEOS-10 (as the pan-Arctic model and analysis/mixed_layer_depth.jl). Buoyancy for comparisons: referenced to the
# surface, b = −g ρ′(Θ, Sᴬ, 0) / ρ₀ (within 1 % of the local value over 256 m); N² in nemo_tke.jl is local.
using Oceananigans.BuoyancyFormulations: SeawaterPolynomials
const teos10 = SeawaterPolynomials.TEOS10.TEOS10EquationOfState()
const reference_density = 1020.0    # SeawaterPolynomials' default for TEOS10EquationOfState
const gravity = 9.80665             # Oceananigans' default
seawater_density_anomaly(T, S, Z) = SeawaterPolynomials.ρ′(T, S, Z, teos10)
seawater_buoyancy(T, S) = -gravity * seawater_density_anomaly(T, S, 0.0) / reference_density

# Horizontally averaged profiles from Oceananigans output (an LES from les.jl or a column from columns.jl), or a NEMO TKE
# column (nemo_tke.jl): centres zc and faces zf (ascending), times t, and b, u, v as vectors of profiles. For seawater
# cases b is computed from T and S (seawater_buoyancy), and T and S are returned too.
function column_profiles(file)
    is_nemo = jldopen(f -> haskey(f, "nemo_tke"), file)
    if is_nemo
        return jldopen(file) do f
            (; zc = f["zc"], zf = f["zf"], t = f["t"], b = f["b"], u = f["u"], v = f["v"],
               T = haskey(f, "T") ? f["T"] : nothing, S = haskey(f, "S") ? f["S"] : nothing)
        end
    end
    seawater = jldopen(f -> haskey(f["timeseries"], "S"), file)
    profiles(fts) = [vec(Array(interior(fts[n], 1, 1, :))) for n in eachindex(fts.times)]
    U, V = FieldTimeSeries(file, "u"), FieldTimeSeries(file, "v")
    if seawater
        Tt, St = FieldTimeSeries(file, "T"), FieldTimeSeries(file, "S")
        T, S = profiles(Tt), profiles(St)
        b = [seawater_buoyancy.(T[n], S[n]) for n in eachindex(T)]
        grid, times = Tt.grid, Tt.times
    else
        B = FieldTimeSeries(file, "b")
        b, T, S = profiles(B), nothing, nothing
        grid, times = B.grid, B.times
    end
    return (; zc = collect(znodes(grid, Center())), zf = collect(znodes(grid, Face())), t = times,
              b, u = profiles(U), v = profiles(V), T, S)
end

# The LES's horizontally averaged output is read the same way
les_profiles(file) = column_profiles(file)

# Vertical grids for the columns: a uniform spacing (m), or the pan-Arctic model's grid ("nz50", "nz100": setup.jl's
# ExponentialDiscretization(Nz, −5000, 0; scale = 1250) cut at the case's depth; its bottom face moved to the domain
# bottom). Faces ascending.
function column_faces(spacing, depth)
    spacing isa Number && return collect(range(-depth, 0; length = round(Int, depth / spacing) + 1))
    Nz = parse(Int, replace(String(spacing), "nz" => ""))
    model_grid = RectilinearGrid(CPU(); size = Nz, z = ExponentialDiscretization(Nz, -5000, 0; scale = 1250),
                                 topology = (Flat, Flat, Bounded))
    faces = collect(znodes(model_grid, Face()))
    faces = faces[faces .> -depth + 1e-6]
    # no sliver at the bottom: the lowest cell at least half the one above it
    length(faces) > 1 && faces[1] + depth < 0.5 * (faces[2] - faces[1]) && (faces = faces[2:end])
    return vcat(-depth, faces)
end
spacing_label(spacing) = spacing isa Number ? "dz$(round(Int, spacing))" : String(spacing)
column_spacings = (1, 4, 8, "nz50", "nz100")
