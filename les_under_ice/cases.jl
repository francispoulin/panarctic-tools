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
#   ice_drag, uᵢ, vᵢ, Cd     drifting ice: instead of τx, the surface stress is quadratic drag on the ice–ocean velocity
#                       difference, τ/ρ₀ = Cd |uᵢ − u₁| (uᵢ − u₁), with the ocean's top-cell velocity u₁ (decision 1 of
#                       the plan: as NEMO SI3 and ClimaSeaIce; Cd = 5e-3)
function buoyancy_case(; τx = 0.0, Jb = 0.0, f = 1e-4, stop_time = 24hours, stokes = false,
                         b_init = three_layer_buoyancy, N²_bottom = three_layer.N²_deep, depth = 256.0,
                         les_spacing = 2.0, les_width = 512.0, ice_drag = false, uᵢ = 0.0, vᵢ = 0.0, Cd = 0.0)
    return (; τx, Jb, Js = 0.0, Jt = 0.0, f, stop_time, stokes, seawater = false,
              b_init, T_init = z -> 0.0, S_init = z -> 0.0, N²_bottom, dTdz_bottom = 0.0, dSdz_bottom = 0.0,
              depth, les_spacing, les_width, ice_drag, uᵢ, vᵢ, Cd)
end

function seawater_case(; τx = 0.0, Js = 0.0, Jt = 0.0, f = 1.4e-4, stop_time = 4days, stokes = false,
                         T_init, S_init, dTdz_bottom, dSdz_bottom, depth = 128.0,
                         les_spacing = 1.0, les_width = 256.0, ice_drag = false, uᵢ = 0.0, vᵢ = 0.0, Cd = 0.0)
    return (; τx, Jb = 0.0, Js, Jt, f, stop_time, stokes, seawater = true,
              b_init = z -> 0.0, T_init, S_init, N²_bottom = 0.0, dTdz_bottom, dSdz_bottom,
              depth, les_spacing, les_width, ice_drag, uᵢ, vᵢ, Cd)
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

# Rung 3: stress without waves (the ice–ocean stress is shear without waves), 24 h, 128 m deep, LES 256² at 1 m.
#   kato_phillips      constant stress u★ = 0.02 m/s (τ = 4e-4 m² s⁻²) into uniform N² = 1e-4 s⁻², no rotation: the
#                      laboratory case with a published deepening law, h = 1.05 u★ (t/N)^½ (Price 1979; used for NEMO's
#                      1-D tests by Reffray et al. 2015): ≈ 62 m at 24 h
#   rotating_stress    the same with f = 1.4e-4 s⁻¹ (75° N): deepening stops after about half an inertial period
#                      (Pollard, Rhines and Thompson 1973)
stress_N² = 1e-4
cases["kato_phillips"]   = buoyancy_case(; τx = -4e-4, f = 0.0, b_init = z -> stress_N² * z, N²_bottom = stress_N²,
                                           depth = 128.0, les_spacing = 1.0, les_width = 256.0)
cases["rotating_stress"] = buoyancy_case(; τx = -4e-4, f = 1.4e-4, b_init = z -> stress_N² * z, N²_bottom = stress_N²,
                                           depth = 128.0, les_spacing = 1.0, les_width = 256.0)

# Rungs 3, 4 and 2 in the Arctic: ice drifting at 0.2 m/s (a storm) over the rung 1 Arctic profile, quadratic drag
# Cd = 5e-3 (initial u★ ≈ 0.014 m/s), 4 days, f = 1.4e-4:
#   arctic_ice_drag        stress only (rung 3)
#   arctic_ice_drag_brine  with the brine of 5 cm/day ice growth, as arctic_brine_convection (rung 4: freeze-up under
#                          drifting ice, the regime of the pan-Arctic bias)
#   arctic_ice_drag_melt   with the freshwater of 2 cm/day melt, Js = +2.3e-7 m/s × 25 g/kg × 0.9 ≈ +5.2e-6 g/kg m/s
#                          (rung 2: a stabilising flux under drifting ice)
arctic_drag = (; ice_drag = true, uᵢ = 0.2, vᵢ = 0.0, Cd = 5e-3)
arctic_initial = (; T_init = z -> arctic_profile.T, S_init = arctic_salinity, dTdz_bottom = 0.0,
                    dSdz_bottom = -arctic_profile.dSdz_halocline)
cases["arctic_ice_drag"]       = seawater_case(; arctic_initial..., arctic_drag...)
cases["arctic_ice_drag_brine"] = seawater_case(; arctic_initial..., arctic_drag..., Js = -1.3e-5)
cases["arctic_ice_drag_melt"]  = seawater_case(; arctic_initial..., arctic_drag..., Js = 5.2e-6)
# Rung 2 proper: melt under calm conditions (ice at rest, no stress, no convection). There is no turbulence source, so
# the truth needs no LES: the freshwater stays where it enters, spreading only by molecular diffusion (√(κ t) ≈ 2 cm
# in 4 days). The test is whether a scheme mixes when it should not (TKE floors, NEMO's surface minimum e ≥ 1e-4,
# background diffusivities). Columns only: calm_melt.jl.
cases["arctic_calm_melt"] = seawater_case(; arctic_initial..., Js = 5.2e-6)
rung3_cases = ["kato_phillips", "rotating_stress", "arctic_ice_drag", "arctic_ice_drag_brine", "arctic_ice_drag_melt"]

# Calibration matrix (23 Nov): 21 LES at 1 m, 4 days (rotating stress 1 day), each one setting away from the cases above,
# so the fit has to follow the scaling with drift speed, latitude, surface flux and stratification. Run by les_matrix.jl.
#   drift speed    uᵢ = 0.1, 0.3 (drag, brine, melt) and 0.4 m/s (drag, brine); 0.2 m/s are the cases above
#   latitude       f = 0 (drag), 0.7e-4 and 1.0e-4 s⁻¹ (drag, brine): is the rotation limit about rotation?
#   halocline      half and double the salinity gradient (drag, brine)
#   mixed layer    20 m instead of 40 m (drag, melt)
#   stress         rotating_stress at f = 0.5e-4 and 1.0e-4 s⁻¹ (with kato_phillips, f = 0, and rotating_stress, 1.4e-4)
# Same initial profile as the cases above (idealised Beaufort winter: 40 m mixed layer at −1.6 °C, S = 30, halocline).
arctic_flux(flux) = flux == :brine ? -1.3e-5 : flux == :melt ? 5.2e-6 : 0.0
function arctic_matrix_case(; uᵢ = 0.2, flux = :drag, f = 1.4e-4, dSdz_halocline = arctic_profile.dSdz_halocline,
                              mixed_layer_depth = arctic_profile.mixed_layer_depth)
    p = merge(arctic_profile, (; dSdz_halocline, mixed_layer_depth))
    return seawater_case(; T_init = z -> p.T, S_init = z -> arctic_salinity(z; p), dTdz_bottom = 0.0,
                           dSdz_bottom = -p.dSdz_halocline, ice_drag = true, uᵢ, vᵢ = 0.0, Cd = 5e-3,
                           Js = arctic_flux(flux), f)
end
flux_suffix(flux) = flux == :drag ? "" : "_" * String(flux)
les_matrix = String[]
function add_matrix_case!(name, case)
    cases[name] = case
    push!(les_matrix, name)
end
for (uᵢ, fluxes) in ((0.1, (:drag, :brine, :melt)), (0.3, (:drag, :brine, :melt)), (0.4, (:drag, :brine)))
    for flux in fluxes
        add_matrix_case!("arctic_ice_drag_u$(uᵢ)" * flux_suffix(flux), arctic_matrix_case(; uᵢ, flux))
    end
end
add_matrix_case!("arctic_ice_drag_f0", arctic_matrix_case(; f = 0.0))
for (label, f) in (("0.7e-4", 0.7e-4), ("1.0e-4", 1.0e-4)), flux in (:drag, :brine)
    add_matrix_case!("arctic_ice_drag_f$(label)" * flux_suffix(flux), arctic_matrix_case(; f, flux))
end
for (label, factor) in (("weak", 0.5), ("strong", 2.0)), flux in (:drag, :brine)
    add_matrix_case!("arctic_ice_drag_$(label)_halocline" * flux_suffix(flux),
                     arctic_matrix_case(; flux, dSdz_halocline = factor * arctic_profile.dSdz_halocline))
end
for flux in (:drag, :melt)
    add_matrix_case!("arctic_ice_drag_ml20" * flux_suffix(flux), arctic_matrix_case(; flux, mixed_layer_depth = 20.0))
end
for (label, f) in (("0.5e-4", 0.5e-4), ("1.0e-4", 1.0e-4))
    add_matrix_case!("rotating_stress_f$(label)",
                     buoyancy_case(; τx = -4e-4, f, b_init = z -> stress_N² * z, N²_bottom = stress_N²,
                                     depth = 128.0, les_spacing = 1.0, les_width = 256.0))
end

# Quadratic ice–ocean drag as top momentum fluxes (positive upward: the ocean gains momentum towards the ice velocity
# when the flux is negative). Discrete form: the top cell's velocity, the other component interpolated to the point.
# LES (3-D, staggered) and column (one point: no interpolation) versions.
using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
@inline function ice_drag_u(i, j, grid, clock, fields, p)
    k = grid.Nz
    Δu = p.uᵢ - @inbounds(fields.u[i, j, k])
    Δv = p.vᵢ - ℑxyᶠᶜᵃ(i, j, k, grid, fields.v)
    return -p.Cd * sqrt(Δu^2 + Δv^2) * Δu
end
@inline function ice_drag_v(i, j, grid, clock, fields, p)
    k = grid.Nz
    Δu = p.uᵢ - ℑxyᶜᶠᵃ(i, j, k, grid, fields.u)
    Δv = p.vᵢ - @inbounds(fields.v[i, j, k])
    return -p.Cd * sqrt(Δu^2 + Δv^2) * Δv
end
@inline function column_ice_drag_u(i, j, grid, clock, fields, p)
    k = grid.Nz
    Δu, Δv = p.uᵢ - @inbounds(fields.u[i, j, k]), p.vᵢ - @inbounds(fields.v[i, j, k])
    return -p.Cd * sqrt(Δu^2 + Δv^2) * Δu
end
@inline function column_ice_drag_v(i, j, grid, clock, fields, p)
    k = grid.Nz
    Δu, Δv = p.uᵢ - @inbounds(fields.u[i, j, k]), p.vᵢ - @inbounds(fields.v[i, j, k])
    return -p.Cd * sqrt(Δu^2 + Δv^2) * Δv
end

# Momentum boundary conditions of a case: constant τx, or ice drag (column = true for single columns)
function momentum_boundary_conditions(case; column = false)
    if case.ice_drag
        p = (; uᵢ = case.uᵢ, vᵢ = case.vᵢ, Cd = case.Cd)
        fu, fv = column ? (column_ice_drag_u, column_ice_drag_v) : (ice_drag_u, ice_drag_v)
        return (u = FieldBoundaryConditions(top = FluxBoundaryCondition(fu, discrete_form = true, parameters = p)),
                v = FieldBoundaryConditions(top = FluxBoundaryCondition(fv, discrete_form = true, parameters = p)))
    else
        return (; u = FieldBoundaryConditions(top = FluxBoundaryCondition(case.τx)))
    end
end

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
               T = haskey(f, "T") ? f["T"] : nothing, S = haskey(f, "S") ? f["S"] : nothing,
               κc = haskey(f, "κc") ? f["κc"] : nothing, κu = haskey(f, "κu") ? f["κu"] : nothing)
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
    # Diffusivities on faces, when the file has them (columns since 22 Nov)
    diffusivity(name) = jldopen(f -> haskey(f["timeseries"], name), file) ? profiles(FieldTimeSeries(file, name)) : nothing
    return (; zc = collect(znodes(grid, Center())), zf = collect(znodes(grid, Face())), t = times,
              b, u = profiles(U), v = profiles(V), T, S, κc = diffusivity("κc"), κu = diffusivity("κu"))
end

# The LES's horizontally averaged output is read the same way
les_profiles(file) = column_profiles(file)

# Vertical grids for the columns: a uniform spacing (m), or the pan-Arctic model's grid ("nz50", "nz100": setup.jl's
# ExponentialDiscretization(Nz, −5000, 0; scale = 1250) cut at the case's depth; its bottom face moved to the domain
# bottom). Faces ascending.
# A function, not a const: cases.jl is included again in reused REPLs.
# "glorys50": GLORYS12's 50 levels (NEMO ORCA standard; T-point depths below), faces halfway between them, top face at 0:
# top cell ~1 m, ~2 m at 10 m depth, ~4 m at 25 m (Greg Wagner on #48, 23 Nov: more vertical resolution, as GLORYS).
glorys_centres() = [0.494, 1.541, 2.646, 3.819, 5.078, 6.441, 7.930, 9.573, 11.405, 13.467, 15.810, 18.496, 21.599,
                        25.211, 29.445, 34.434, 40.344, 47.374, 55.764, 65.807, 77.854, 92.326, 109.729, 130.666, 155.851,
                        186.126, 222.475, 266.040, 318.127, 380.213, 453.938, 541.089, 643.567, 763.333, 902.339,
                        1062.440, 1245.291, 1452.251, 1684.284, 1941.893, 2225.078, 2533.336, 2865.703, 3220.820,
                        3597.032, 3992.484, 4405.224, 4833.291, 5274.784, 5727.917]
function column_faces(spacing, depth)
    spacing isa Number && return collect(range(-depth, 0; length = round(Int, depth / spacing) + 1))
    if String(spacing) == "glorys50"
        faces = sort(-vcat(0.0, (glorys_centres()[1:end-1] .+ glorys_centres()[2:end]) ./ 2))
    else
        Nz = parse(Int, replace(String(spacing), "nz" => ""))
        model_grid = RectilinearGrid(CPU(); size = Nz, z = ExponentialDiscretization(Nz, -5000, 0; scale = 1250),
                                     topology = (Flat, Flat, Bounded))
        faces = collect(znodes(model_grid, Face()))
    end
    faces = faces[faces .> -depth + 1e-6]
    # no sliver at the bottom: the lowest cell at least half the one above it
    length(faces) > 1 && faces[1] + depth < 0.5 * (faces[2] - faces[1]) && (faces = faces[2:end])
    return vcat(-depth, faces)
end
spacing_label(spacing) = spacing isa Number ? "dz$(round(Int, spacing))" : String(spacing)
column_spacings = (1, 4, 8, "nz50", "nz100")

# Base of the mixed layer for a uniform background N²: the shallowest depth below 5 m where N² exceeds half the
# background (searching downward; zc ascending). `interpolate = true` places it between faces by linear interpolation
# of N² across the threshold, so coarse grids are not limited to whole cells (rung 1: nz50 cells are 11–14 m at 70–100 m).
function mixed_layer_base(b, zc, depth, N²_background; interpolate = false)
    threshold = 0.5N²_background
    previous_N², previous_z = NaN, NaN
    for k in length(b)-1:-1:1
        zf = (zc[k] + zc[k+1]) / 2
        -zf < 5 && continue
        -zf > 0.8depth && break
        N² = (b[k+1] - b[k]) / (zc[k+1] - zc[k])
        if N² > threshold
            (interpolate && isfinite(previous_N²)) || return -zf
            w = (threshold - previous_N²) / (N² - previous_N²)
            return -(previous_z + w * (zf - previous_z))
        end
        previous_N², previous_z = N², zf
    end
    return NaN
end

# The thickness of the column cell containing depth h (for the "one cell" part of the pass criterion)
cell_thickness_at(h, zf) = (k = findlast(z -> z ≤ -h, zf); k === nothing || k == length(zf) || isnan(h) ? 0.0 : zf[k+1] - zf[k])

# Mixed-layer depth from a buoyancy profile b on centres zc (ascending), as in compare.jl: the first depth below z_ref
# where b is below its value at z_ref by more than Δb. interpolate = true places it between centres (smooth in the
# parameters, for screening and calibration).
function threshold_depth(b, z; Δb = 1e-5, z_ref = -10.0, interpolate = false)
    reference = b[argmin(abs.(z .- z_ref))]
    k = findlast(k -> z[k] < z_ref && b[k] < reference - Δb, eachindex(b))
    k === nothing && return -z[1]
    (interpolate && k < length(b) && b[k+1] != b[k]) || return -z[k]
    w = (reference - Δb - b[k+1]) / (b[k] - b[k+1])
    return -(z[k+1] + w * (z[k] - z[k+1]))
end

# The depth measure of each kind of case: the base of the mixed layer for uniform N² (stress and convection cases), the
# threshold depth otherwise
uniform_N²(name) = (name == "kato_phillips" || startswith(name, "rotating_stress")) ? stress_N² :
                   startswith(name, "convection_") ? first(p[1] for p in convection_suite if convection_case_name(p...) == name) :
                   nothing
case_depth(name, case, b, zc; interpolate = true) =
    (N² = uniform_N²(name);
     N² === nothing ? threshold_depth(b, zc; interpolate) : mixed_layer_base(b, zc, case.depth, N²; interpolate))
