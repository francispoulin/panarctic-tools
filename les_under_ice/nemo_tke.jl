# NEMO's TKE vertical mixing scheme (zdftke.F90, NEMO 4.x; Gaspar et al. 1990, Blanke and Delecluse 1993, Madec et al.)
# in a standalone single column, with the same cases, forcing and output as columns.jl, so compare.jl and
# compare_lesbrary.jl show it next to CATKE. It separates what CATKE folds together (claude/catke-les-plan.md):
#   base scheme   TKE equation with shear and buoyancy production, dissipation ediss e^{3/2}/ℓ, diffusion; surface
#                 value e = ebb |τ| (also under ice, from the ice–ocean stress); mixing length √(2e)/N bounded by the
#                 distance to the boundaries (nn_mxl = 2), surface length from the stress (ln_mxl0); Prandtl number
#                 from the Richardson number (nn_pdl = 1)
#   waves         Langmuir cells (ln_lc: a TKE source over a depth set by the potential energy of the stratification
#                 against the Langmuir velocity, parameterised as 0.016 U₁₀) and penetration of surface TKE below the
#                 mixed layer (nn_etau = 1: rn_efr e_surface exp(−z/htau)); both × (1 − w(ice fraction)), with
#                 w = tanh(10 fᵢ) by default (nn_eice = 1, our "tanh")
# Defaults are NEMO's reference namelist (namzdf_tke): rn_ediff = 0.1, rn_ediss = 0.7, rn_ebb = 67.83, rn_emin = 1e-6,
# rn_emin0 = 1e-4, nn_mxl = 2, ln_mxl0 = true, nn_pdl = 1, ln_lc = true, rn_lc = 0.15, nn_etau = 1, rn_efr = 0.05,
# nn_htau = 1, nn_eice = 1; background viscosity and diffusivity 1.2e-4 and 1.2e-5 m²/s, and enhanced vertical
# diffusion of tracers where N² < 0 (namzdf: ln_zdfevd = true, nn_evdm = 0, rn_evd = 100 m²/s). NEMO convects through
# this, not through the TKE: without it (checked in a numpy port of this file, 21 Nov) free convection leaves the
# top 70 m statically unstable after 24 h; with it the convective layer is 102 m at 4 m (our LES: 107 m).
#
# Variants written by this script, for a case (CPU, seconds each):
#   nemo          the reference namelist (Langmuir cells and penetration on)
#   nemo_no_waves base scheme only (ln_lc = false, nn_etau = 0): NEMO under compact ice, or wind without waves
#   nemo_lc_only  Langmuir cells on, penetration off
# each at 1, 4, 8 m and on the pan-Arctic grids at 50 and 100 levels; buoyancy cases or seawater (T, S, TEOS-10) cases.
#
#     empty!(ARGS); push!(ARGS, "strong_wind")
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/nemo_tke.jl")
#
# Internally top-down (index 1 = surface), as NEMO: cells j = 1…N, faces k = 1…N+1 (TKE and diffusivities on faces).
# Each step, in NEMO's order: TKE (Langmuir source added, implicit solve with semi-implicit dissipation, penetration),
# mixing length and diffusivities from the new TKE, then u, v, b by implicit vertical diffusion with the surface
# fluxes, and an exact Coriolis rotation. Written from zdftke.F90 and namzdf_tke; not checked against NEMO's own 1-D
# configuration (C1D).

include(joinpath(@__DIR__, "cases.jl"))
using JLD2, Printf

Base.@kwdef struct NEMOTKEParameters
    rn_ediff  :: Float64 = 0.1
    rn_ediss  :: Float64 = 0.7
    rn_ebb    :: Float64 = 67.83
    rn_emin   :: Float64 = 1e-6
    rn_emin0  :: Float64 = 1e-4
    rn_bshear :: Float64 = 1e-20
    ln_mxl0   :: Bool    = true
    ln_lc     :: Bool    = true
    rn_lc     :: Float64 = 0.15
    nn_etau   :: Int     = 1
    rn_efr    :: Float64 = 0.05
    nn_eice   :: Int     = 1
    ln_zdfevd :: Bool    = true
    nn_evdm   :: Int     = 0
    rn_evd    :: Float64 = 100.0
    ice_fraction :: Float64 = 0.0
    avm0      :: Float64 = 1.2e-4
    avt0      :: Float64 = 1.2e-5
    rho0      :: Float64 = 1026.0
    grav      :: Float64 = 9.80665
    vkarmn    :: Float64 = 0.4
    rho_air   :: Float64 = 1.22
    Cdrag     :: Float64 = 1.5e-3
end

function ice_weight(p::NEMOTKEParameters)
    fᵢ = p.ice_fraction
    p.nn_eice == 0 && return 0.0
    p.nn_eice == 1 && return tanh(10fᵢ)
    p.nn_eice == 2 && return fᵢ
    p.nn_eice == 3 && return min(4fᵢ, 1.0)
    error("nn_eice must be 0, 1, 2 or 3")
end

# Thomas algorithm: a (sub), d (diagonal), c (super), r (right-hand side); a[1] and c[end] unused
function tridiagonal_solve(a, d, c, r)
    n = length(d)
    c′, r′ = zeros(n), zeros(n)
    c′[1], r′[1] = c[1] / d[1], r[1] / d[1]
    for i in 2:n
        m = d[i] - a[i] * c′[i-1]
        c′[i] = i < n ? c[i] / m : 0.0
        r′[i] = (r[i] - a[i] * r′[i-1]) / m
    end
    x = zeros(n)
    x[n] = r′[n]
    for i in n-1:-1:1
        x[i] = r′[i] - c′[i] * x[i+1]
    end
    return x
end

# Implicit vertical diffusion of a cell variable with diffusivities κ on faces (top-down: face j above cell j), cell
# thicknesses Δc, distances between centres Δw (on faces), and a surface tendency
function diffuse(ϕ, κ, Δc, Δw, Δt, surface_tendency)
    N = length(ϕ)
    a, d, c = zeros(N), ones(N), zeros(N)
    for j in 1:N
        j > 1 && (a[j] = -Δt * κ[j] / (Δc[j] * Δw[j]); d[j] -= a[j])
        j < N && (c[j] = -Δt * κ[j+1] / (Δc[j] * Δw[j+1]); d[j] -= c[j])
    end
    r = copy(ϕ)
    r[1] += Δt * surface_tendency
    return tridiagonal_solve(a, d, c, r)
end

# N² on the interior faces k = 2…N, between cell k−1 above and cell k below: from b, or for seawater from the
# density difference at the face's depth (TEOS-10, locally referenced)
function stratification(state, dw, Δw, seawater)
    N = length(dw) - 1
    N² = zeros(N + 1)
    for k in 2:N
        if seawater
            Z = -dw[k]
            Δρ = seawater_density_anomaly(state.T[k], state.S[k], Z) - seawater_density_anomaly(state.T[k-1], state.S[k-1], Z)
            N²[k] = gravity * Δρ / (reference_density * Δw[k])
        else
            N²[k] = (state.b[k-1] - state.b[k]) / Δw[k]
        end
    end
    return N²
end

# A column on the faces of cases.jl's column_faces (uniform spacing in m, or "nz50" / "nz100": the pan-Arctic grids)
function nemo_tke_column(case_name, case, p::NEMOTKEParameters; spacing = 4, Δt = 60.0, label = "nemo",
                         output_interval = 600.0)
    dw = reverse(-column_faces(spacing, case.depth))    # face depths (positive down), top-down; dw[1] = 0
    N  = length(dw) - 1
    dt = [(dw[j] + dw[j+1]) / 2 for j in 1:N]           # cell-centre depths
    Δc = diff(dw)                                       # cell thicknesses (e3t)
    Δw = vcat(dt[1], diff(dt), dw[N+1] - dt[N])         # distances between centres, on faces (e3w)
    seawater = case.seawater
    b = seawater ? zeros(N) : case.b_init.(-dt)
    T = seawater ? case.T_init.(-dt) : zeros(N)
    S = seawater ? case.S_init.(-dt) : zeros(N)
    u, v = zeros(N), zeros(N)
    en = fill(p.rn_emin, N + 1)
    avm, avt = fill(p.avm0, N + 1), fill(p.avt0, N + 1)
    avm[1] = avm[N+1] = avt[1] = avt[N+1] = 0.0

    zice = ice_weight(p)
    ri_cri   = 2 / (2 + p.rn_ediss / p.rn_ediff)
    rmxl_min = 1e-6 / (p.rn_ediff * sqrt(p.rn_emin))
    rn_mxl0  = p.ln_mxl0 ? rmxl_min : 0.04
    zraug    = p.vkarmn * 2e5 / (p.rho0 * p.grav)
    sinφ     = case.f / (2 * 7.292115e-5)
    htau     = max(0.5, min(30.0, 45 * abs(sinφ)))     # nn_htau = 1

    t = 0.0
    times, bs, Ts, Ss, us, vs = [t], [copy(b)], [copy(T)], [copy(S)], [copy(u)], [copy(v)]
    seawater && (bs[1] = seawater_buoyancy.(T, S))
    next_output = output_interval
    nsteps = round(Int, case.stop_time / Δt)
    for _ in 1:nsteps
        # Surface stress (kinematic, positive upward as an Oceananigans flux): constant, or quadratic drag on drifting
        # ice from the top cell's velocity (cases.jl); |τ| in N/m² for the surface TKE, length and Langmuir terms
        if case.ice_drag
            Δu, Δv = case.uᵢ - u[1], case.vᵢ - v[1]
            τx, τy = -case.Cd * hypot(Δu, Δv) * Δu, -case.Cd * hypot(Δu, Δv) * Δv
        else
            τx, τy = case.τx, 0.0
        end
        taum  = p.rho0 * hypot(τx, τy)
        zWlc2 = 0.5 * 0.016^2 / (p.rho_air * p.Cdrag) * taum   # (Langmuir velocity)² / 2, uncoupled
        zus3  = (1 - zice) * sqrt(2zWlc2)^3

        # Stratification and shear on interior faces k = 2…N (between cells k−1 above and k below)
        N² = stratification((; b, T, S), dw, Δw, seawater)
        S² = zeros(N + 1)
        for k in 2:N
            S²[k] = ((u[k-1] - u[k])^2 + (v[k-1] - v[k])^2) / Δw[k]^2
        end
        sh2 = avm .* S²

        # Mixing length (dissipation) from the current TKE, for the dissipation term
        ℓ = mixing_length(en, N², taum, p, rmxl_min, rn_mxl0, zraug, Δc)
        dissl = sqrt.(en) ./ ℓ

        # Surface and bottom values
        en[1]   = max(p.rn_emin0, p.rn_ebb / p.rho0 * taum)
        en[N+1] = p.rn_emin

        # Langmuir cells: depth where the potential energy of the stratification exceeds zWlc2
        if p.ln_lc && zWlc2 > 0
            zpelc = 0.0
            imlc = N + 1
            pe = zeros(N + 1)
            for k in 2:N+1
                zpelc += max(N²[k], 0.0) * dw[k] * Δw[k]
                pe[k] = zpelc
            end
            for k in N:-1:2
                pe[k] > zWlc2 && (imlc = k)                    # the shallowest such face
            end
            zhlc = dw[imlc]
            for k in 2:N
                if dw[k] < zhlc
                    zwlc = p.rn_lc * sin(π * dw[k] / zhlc)
                    en[k] += Δt * zus3 * zwlc^3 / zhlc
                end
            end
        end

        # Implicit TKE equation on interior faces k = 2…N (NEMO: e3t(k) is the cell below face k, here Δc[k])
        n = N - 1
        a, d, c, r = zeros(n), ones(n), zeros(n), zeros(n)
        for (i, k) in enumerate(2:N)
            zd_up = -0.5Δt * max(avm[k+1] + avm[k], 2e-5) / (Δc[k] * Δw[k])
            zd_lw = -0.5Δt * max(avm[k] + avm[k-1], 2e-5) / (Δc[k-1] * Δw[k])
            a[i], c[i] = zd_lw, zd_up
            d[i] = 1 - zd_lw - zd_up + 1.5Δt * p.rn_ediss * dissl[k]
            r[i] = en[k] + Δt * (sh2[k] - avt[k] * N²[k] + 0.5 * p.rn_ediss * dissl[k] * en[k])
        end
        r[1] -= a[1] * en[1]
        r[n] -= c[n] * en[N+1]
        en[2:N] .= tridiagonal_solve(a, d, c, r)

        # Penetration of surface TKE below the mixed layer
        if p.nn_etau == 1
            for k in 2:N
                en[k] += p.rn_efr * en[1] * exp(-dw[k] / htau) * (1 - zice)
            end
        end
        en .= max.(en, p.rn_emin)

        # Diffusivities from the new TKE (tke_avn)
        ℓ = mixing_length(en, N², taum, p, rmxl_min, rn_mxl0, zraug, Δc)
        for k in 2:N
            zav = p.rn_ediff * ℓ[k] * sqrt(en[k])
            avm[k] = max(zav, p.avm0)
            avt[k] = max(zav, p.avt0)
            zri = N²[k] > 0 ? N²[k] * avm[k] / (sh2[k] + p.rn_bshear) : 0.0   # as NEMO: sh2 from the previous avm
            apdlr = max(0.1, ri_cri / max(ri_cri, zri))
            avt[k] = max(apdlr * avt[k], p.avt0)
            if p.ln_zdfevd && N²[k] < 0                       # zdf_evd: convection
                avt[k] = p.rn_evd
                p.nn_evdm == 1 && (avm[k] = p.rn_evd)
            end
        end

        # Momentum and tracers: implicit vertical diffusion with surface fluxes (flux positive upward)
        u = diffuse(u, avm, Δc, Δw, Δt, -τx / Δc[1])
        v = diffuse(v, avm, Δc, Δw, Δt, -τy / Δc[1])
        if seawater
            T = diffuse(T, avt, Δc, Δw, Δt, -case.Jt / Δc[1])
            S = diffuse(S, avt, Δc, Δw, Δt, -case.Js / Δc[1])
        else
            b = diffuse(b, avt, Δc, Δw, Δt, -case.Jb / Δc[1])
        end
        θ = case.f * Δt
        u, v = u .* cos(θ) .+ v .* sin(θ), -u .* sin(θ) .+ v .* cos(θ)

        t += Δt
        if t ≥ next_output - 1e-6
            push!(times, t); push!(us, copy(u)); push!(vs, copy(v)); push!(Ts, copy(T)); push!(Ss, copy(S))
            push!(bs, seawater ? seawater_buoyancy.(T, S) : copy(b))
            next_output += output_interval
        end
    end

    # Saved bottom-up (ascending z), as Oceananigans
    file = joinpath(output_folder(), @sprintf("%s_column_%s_%s.jld2", case_name, label, spacing_label(spacing)))
    jldopen(file, "w") do f
        f["nemo_tke"] = true
        f["parameters"] = Dict(String(n) => getfield(p, n) for n in fieldnames(NEMOTKEParameters))
        f["zc"] = reverse(-dt)
        f["zf"] = reverse(-dw)
        f["t"] = times
        f["b"] = reverse.(bs); f["u"] = reverse.(us); f["v"] = reverse.(vs)
        if seawater
            f["T"] = reverse.(Ts); f["S"] = reverse.(Ss)
        end
    end
    @info "NEMO TKE column $case_name, $label, $(spacing_label(spacing)): $file"
    return file
end

# nn_mxl = 2: ℓ = √(2e)/N, surface value from the stress (ln_mxl0), limited by the distance to the boundaries (growing
# by at most the cell thickness from one face to the next)
function mixing_length(en, N², taum, p, rmxl_min, rn_mxl0, zraug, Δc)
    K = length(en)
    ℓ = [max(rmxl_min, sqrt(2en[k] / max(N²[k], 1e-20))) for k in 1:K]
    ℓ[1] = p.ln_mxl0 ? max(rn_mxl0, zraug * taum) : rn_mxl0
    ℓ[K] = rmxl_min
    for k in 2:K-1
        ℓ[k] = min(ℓ[k-1] + Δc[k-1], ℓ[k])
    end
    for k in K-1:-1:2
        ℓ[k] = min(ℓ[k+1] + Δc[k], ℓ[k])
    end
    return ℓ
end

nemo_case = isempty(ARGS) ? "strong_wind" : ARGS[1]
# nemo_no_breaking (rung 3, 22 Nov): no wave terms and no wave breaking either: Gaspar et al.'s surface value
# e = 3.75 u★² instead of 67.83 u★², and the fixed 0.04 m surface length instead of the Charnock-type one. NEMO under
# ice keeps 67.83 (its ice weighting does not touch it), so this asks whether that matters under drifting ice.
nemo_variants = Dict("nemo"             => NEMOTKEParameters(),
                     "nemo_no_waves"    => NEMOTKEParameters(ln_lc = false, nn_etau = 0),
                     "nemo_lc_only"     => NEMOTKEParameters(nn_etau = 0),
                     "nemo_no_breaking" => NEMOTKEParameters(ln_lc = false, nn_etau = 0, rn_ebb = 3.75, ln_mxl0 = false))
for (label, parameters) in sort(collect(nemo_variants); by = first), spacing in column_spacings
    nemo_tke_column(nemo_case, case_parameters(nemo_case), parameters; spacing, label)
end
