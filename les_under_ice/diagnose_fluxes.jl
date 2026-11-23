# Step 1 of the CATKE recalibration plan (claude/catke-les-plan.md, 22 Nov): where do the schemes over-mix?
#
# The LES saved horizontally averaged turbulent fluxes (⟨w b⟩ or ⟨w T⟩, ⟨w S⟩; ⟨w u⟩, ⟨w v⟩). Divided by the mean gradients
# they give the LES's own eddy diffusivities, as functions of depth:
#   κc = −⟨w b⟩ / ∂z b                                   (tracers; for seawater ⟨w b⟩ = −(g/ρ₀)(ρ_T ⟨w T⟩ + ρ_S ⟨w S⟩))
#   κu = −(⟨w u⟩ ∂z u + ⟨w v⟩ ∂z v) / |∂z u|²            (momentum, the flux projected on the shear)
#   Ri = N² / |∂z u|²,  Pr = κu / κc
# compared with the diffusivities of the columns (CATKE and k-ε at 1 and 4 m, NEMO base and NEMO without breaking at
# 1 m), which are rerun here so that they save κu as well as κc (CPU, a few minutes per case).
#
# Each quantity is computed at every output time in the second half of the run and its median is taken (robust to
# inertial oscillations and noise); κ is left out where the gradient is weak (N² or S² below 2 % of the profile's
# maximum). The LES fluxes are resolved fluxes only (implicit LES): too small within a few metres of the surface.
#
# Writes output/diagnose_<case>.png and output/diagnose_fluxes.md: κc and κu in two bands relative to the LES mixed-layer
# depth h(t) (interior 0.3–0.8 h, base 0.8–1.1 h; 23 Nov: the first version sampled one depth below the layer, in
# quiet water), for the LES and each column, and the ratios to the LES.
#
#     empty!(ARGS)                       # or case names; default: kato_phillips, rotating_stress, arctic_ice_drag,
#                                        # med_wind_med_cooling
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/diagnose_fluxes.jl")

using Oceananigans, Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity, TKEDissipationVerticalDiffusivity
using Oceananigans.OutputReaders: FieldTimeSeries
using CairoMakie, JLD2, Printf, Statistics

include(joinpath(@__DIR__, "cases.jl"))
include(joinpath(@__DIR__, "column_model.jl"))
include(joinpath(@__DIR__, "nemo_model.jl"))

diagnosed_cases = isempty(ARGS) ? ["kato_phillips", "rotating_stress", "arctic_ice_drag", "med_wind_med_cooling"] : copy(ARGS)
empty!(ARGS)

# The columns compared (name in the figure => how to (re)run it)
diagnosed_columns = [("CATKE 1 m", "stock", 1), ("CATKE 4 m", "stock", 4), ("k-ε 1 m", "k_epsilon", 1),
                     ("k-ε 4 m", "k_epsilon", 4), ("NEMO 1 m", "nemo_no_waves", 1), ("NEMO no breaking 1 m", "nemo_no_breaking", 1)]
column_closures = Dict("stock" => CATKEVerticalDiffusivity(), "k_epsilon" => TKEDissipationVerticalDiffusivity())

nanmedian(x) = (y = filter(isfinite, x); isempty(y) ? NaN : median(y))

# A profile of fluxes onto the interior faces (between centres): from centres (average) or faces (drop the ends)
function to_interior_faces(values, Nz)
    length(values) == Nz     && return (values[1:end-1] .+ values[2:end]) ./ 2
    length(values) == Nz + 1 && return values[2:end-1]
    error("Flux profile of length $(length(values)) for Nz = $Nz")
end

profile_series(file, name) = (fts = FieldTimeSeries(file, name); [vec(Array(interior(fts[n], 1, 1, :))) for n in eachindex(fts.times)])

# ρ_T and ρ_S at the surface (as seawater_buoyancy), by centred differences
density_derivatives(T, S) = ((seawater_density_anomaly(T + 1e-3, S, 0.0) - seawater_density_anomaly(T - 1e-3, S, 0.0)) / 2e-3,
                             (seawater_density_anomaly(T, S + 1e-3, 0.0) - seawater_density_anomaly(T, S - 1e-3, 0.0)) / 2e-3)

# Profiles on interior faces at each output time of the second half of the run: Ri, κc, κu (NaN where the gradient is
# weak: N² or S² below 2 % of the profile's maximum). Returned as matrices (face × time) with the times.
function diffusivity_matrices(zc, t, b, u, v, flux_b, flux_u, flux_v)
    zf = (zc[1:end-1] .+ zc[2:end]) ./ 2
    Δz = diff(zc)
    late = findall(≥(t[end] / 2), t)
    Ri, κc, κu = (fill(NaN, length(zf), length(late)) for _ in 1:3)
    for (m, n) in enumerate(late)
        N² = diff(b[n]) ./ Δz
        S² = (diff(u[n]) ./ Δz) .^ 2 .+ (diff(v[n]) ./ Δz) .^ 2
        strong_N², strong_S² = N² .> 0.02maximum(N²), S² .> 0.02maximum(S²)
        Ri[:, m] = ifelse.(strong_S², N² ./ S², NaN)
        if flux_b === nothing                    # a column: its own diffusivities
            κc[:, m] = flux_u[n]
            κu[:, m] = flux_v[n]
        else
            κc[:, m] = ifelse.(strong_N², -flux_b[n] ./ N², NaN)
            momentum = -(flux_u[n] .* diff(u[n]) ./ Δz .+ flux_v[n] .* diff(v[n]) ./ Δz)
            κu[:, m] = ifelse.(strong_S², momentum ./ S², NaN)
        end
    end
    return (; zf, t = t[late], Ri, κc, κu)
end

function les_diffusivities(file, case)
    les = les_profiles(file)
    Nz = length(les.zc)
    wu, wv = profile_series(file, "wu"), profile_series(file, "wv")
    if case.seawater
        wT, wS = profile_series(file, "wT"), profile_series(file, "wS")
        wb = map(eachindex(les.t)) do n
            T, S = to_interior_faces(les.T[n], Nz), to_interior_faces(les.S[n], Nz)
            derivatives = density_derivatives.(T, S)
            ρT, ρS = first.(derivatives), last.(derivatives)
            -(gravity / reference_density) .* (ρT .* to_interior_faces(wT[n], Nz) .+ ρS .* to_interior_faces(wS[n], Nz))
        end
    else
        wb = [to_interior_faces(f, Nz) for f in profile_series(file, "wb")]
    end
    fu = [to_interior_faces(f, Nz) for f in wu]
    fv = [to_interior_faces(f, Nz) for f in wv]
    return les, diffusivity_matrices(les.zc, les.t, les.b, les.u, les.v, wb, fu, fv)
end

# A column's own diffusivities on its interior faces (saved on all faces), as matrices like the LES's
column_diffusivities(col) =
    diffusivity_matrices(col.zc, col.t, col.b, col.u, col.v, nothing,
                         [κ[2:end-1] for κ in col.κc],
                         col.κu === nothing ? [fill(NaN, length(col.zc) - 1) for _ in col.κc] : [κ[2:end-1] for κ in col.κu])

# Median of a quantity over the band lo h(t) < depth < hi h(t), all times of the second half; h(t) from the LES
function band_median(matrix, zf, times, h_of_time, lo, hi)
    values = Float64[]
    for (m, t) in enumerate(times)
        h = h_of_time(t)
        for k in eachindex(zf)
            lo * h < -zf[k] < hi * h && push!(values, matrix[k, m])
        end
    end
    return nanmedian(values)
end
row_median(matrix) = map(nanmedian, eachrow(matrix))

# Two bands, relative to the LES mixed-layer depth h(t): the interior of the boundary layer (0.3–0.8 h), where the
# resolved fluxes are reliable, and the base (0.8–1.1 h), where entrainment happens
bands = (interior = (0.3, 0.8), base = (0.8, 1.1))

function diagnose(case_name, folder, rows)
    case = case_parameters(case_name)
    les_file = joinpath(folder, case_name * "_les.jld2")
    isfile(les_file) || (@warn "No LES for $case_name: skipped"; return)
    les, D = les_diffusivities(les_file, case)
    h_les = [case_depth(case_name, case, b, les.zc) for b in les.b]
    h_of_time(t) = h_les[argmin(abs.(les.t .- t))]
    h_end = h_les[end]

    band_values(M) = ([band_median(M.κc, M.zf, M.t, h_of_time, bands[k]...) for k in keys(bands)],
                      [band_median(M.κu, M.zf, M.t, h_of_time, bands[k]...) for k in keys(bands)],
                      band_median(M.Ri, M.zf, M.t, h_of_time, bands.base...))
    κc_les, κu_les, Ri_les = band_values(D)
    push!(rows, @sprintf("| %s | LES | %.1f | %.2f | %.1e | %.1e | %.1e | %.1e | | | | |", case_name, h_end, Ri_les,
                         κc_les[1], κc_les[2], κu_les[1], κu_les[2]))

    fig = Figure(size = (1600, 650))
    titles = ["Ri = N² / S²", "κc (m²/s), tracers", "κu (m²/s), momentum", "Pr = κu / κc"]
    panels = [Axis(fig[1, i]; title = titles[i], ylabel = i == 1 ? "z (m)" : "", xscale = i == 1 ? identity : log10)
              for i in 1:4]
    depth_shown = -min(0.6case.depth, 2h_end)
    plot_profile!(ax, x, z; kw...) = (keep = isfinite.(x) .& (x .> 0) .& (z .> depth_shown);
                                      any(keep) && lines!(ax, x[keep], z[keep]; kw...))
    plot_all!(M; kw...) = (plot_profile!(panels[1], clamp.(row_median(M.Ri), 1e-3, 10), M.zf; kw...);
                           plot_profile!(panels[2], row_median(M.κc), M.zf; kw...);
                           plot_profile!(panels[3], row_median(M.κu), M.zf; kw...);
                           plot_profile!(panels[4], row_median(M.κu) ./ row_median(M.κc), M.zf; kw...))
    plot_all!(D; color = :black, linewidth = 3, label = "LES")

    colours = Makie.wong_colors()
    for (c, (label, variant, spacing)) in enumerate(diagnosed_columns)
        if startswith(variant, "nemo")
            file = nemo_tke_column(case_name, case, nemo_variants[variant]; spacing, label = variant)
        else
            file = run_column(case_name, case, variant, column_closures[variant], spacing)
        end
        col = column_profiles(file)
        C = column_diffusivities(col)
        plot_all!(C; color = colours[c], linestyle = spacing == 1 ? :solid : :dash, label)
        κc_col, κu_col, Ri_col = band_values(C)
        h_col = case_depth(case_name, case, col.b[end], col.zc)
        push!(rows, @sprintf("| %s | %s | %.1f | %.2f | %.1e | %.1e | %.1e | %.1e | %.1f | %.1f | %.1f | %.1f |",
                             case_name, label, h_col, Ri_col, κc_col[1], κc_col[2], κu_col[1], κu_col[2],
                             κc_col[1] / κc_les[1], κc_col[2] / κc_les[2], κu_col[1] / κu_les[1], κu_col[2] / κu_les[2]))
    end
    for ax in panels
        hlines!(ax, [-h_end, -0.3h_end]; color = :gray, linestyle = :dot)
        ylims!(ax, depth_shown, 0)
    end
    axislegend(panels[1]; position = :rb, labelsize = 9)
    Label(fig[0, :], "$case_name: diffusivities from the LES fluxes (black) and the columns, median over the second half " *
                     "(dotted: 0.3 h and h of the LES at the end)"; fontsize = 15)
    save(joinpath(folder, "diagnose_$(case_name).png"), fig)
    return nothing
end

function diagnose_report(case_names, folder)
    rows = String[]
    for name in case_names
        diagnose(name, folder, rows)
    end
    io = IOBuffer()
    println(io, "# Eddy diffusivities in the boundary layer: LES fluxes against the columns\n")
    println(io, "Medians over the second half of the run, in two bands relative to the LES mixed-layer depth h(t): interior ",
                "0.3–0.8 h and base 0.8–1.1 h. κ_LES from the resolved fluxes (negative or tiny values: no resolved mixing). ",
                "h: each profile's own depth at the end. Ri at the base band. Ratios > 1: the column mixes more than the LES.\n")
    println(io, "| Case | Profile | h (m) | Ri base | κc interior | κc base | κu interior | κu base | κc/LES interior | ",
                "κc/LES base | κu/LES interior | κu/LES base |")
    println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|")
    foreach(r -> println(io, r), rows)
    text = String(take!(io))
    print(text)
    write(joinpath(folder, "diagnose_fluxes.md"), text)
    return text
end

diagnose_report(diagnosed_cases, output_folder())
@info "Wrote $(joinpath(output_folder(), "diagnose_fluxes.md")) and diagnose_<case>.png"
