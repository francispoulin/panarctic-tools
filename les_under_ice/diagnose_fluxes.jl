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
# Writes output/diagnose_<case>.png and output/diagnose_fluxes.md: at the base of the mixed layer in the LES (the depth
# of the maximum of N²), Ri, κc, κu and Pr for the LES and each column, and each column's κ / κ_LES.
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

# Median profiles over the second half of the run on interior faces: Ri, κc, κu (NaN where the gradient is weak)
function median_diffusivities(zc, t, b, u, v, flux_b, flux_u, flux_v)
    zf = (zc[1:end-1] .+ zc[2:end]) ./ 2
    Δz = diff(zc)
    late = findall(≥(t[end] / 2), t)
    Ri, κc, κu = (fill(NaN, length(zf), length(late)) for _ in 1:3)
    for (m, n) in enumerate(late)
        N² = diff(b[n]) ./ Δz
        S² = (diff(u[n]) ./ Δz) .^ 2 .+ (diff(v[n]) ./ Δz) .^ 2
        strong_N², strong_S² = N² .> 0.02maximum(N²), S² .> 0.02maximum(S²)
        Ri[:, m] = ifelse.(strong_S², N² ./ S², NaN)
        κc[:, m] = ifelse.(strong_N², -flux_b[n] ./ N², NaN)
        momentum = -(flux_u[n] .* diff(u[n]) ./ Δz .+ flux_v[n] .* diff(v[n]) ./ Δz)
        κu[:, m] = ifelse.(strong_S², momentum ./ S², NaN)
    end
    return (; zf, Ri = map(nanmedian, eachrow(Ri)), κc = map(nanmedian, eachrow(κc)), κu = map(nanmedian, eachrow(κu)))
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
    return median_diffusivities(les.zc, les.t, les.b, les.u, les.v, wb, fu, fv)
end

# Column diffusivities are saved on all faces zf (ascending): keep the interior ones, median over the second half
function column_diffusivities(col)
    zf = col.zf[2:end-1]
    late = findall(≥(col.t[end] / 2), col.t)
    κc = [nanmedian([col.κc[n][k+1] for n in late]) for k in eachindex(zf)]
    κu = col.κu === nothing ? fill(NaN, length(zf)) : [nanmedian([col.κu[n][k+1] for n in late]) for k in eachindex(zf)]
    Δz = diff(col.zc)
    Ri = map(eachindex(zf)) do k
        values = map(late) do n
            N² = (col.b[n][k+1] - col.b[n][k]) / Δz[k]
            S² = ((col.u[n][k+1] - col.u[n][k])^2 + (col.v[n][k+1] - col.v[n][k])^2) / Δz[k]^2
            S² > 0 ? N² / S² : NaN
        end
        nanmedian(values)
    end
    return (; zf, Ri, κc, κu)
end

nearest_value(zf, values, z) = values[argmin(abs.(zf .- z))]

function diagnose(case_name, folder, rows)
    case = case_parameters(case_name)
    les_file = joinpath(folder, case_name * "_les.jld2")
    isfile(les_file) || (@warn "No LES for $case_name: skipped"; return)
    les = les_diffusivities(les_file, case)

    # Base of the mixed layer in the LES: the depth of the maximum of N², from the final profile
    les_end = les_profiles(les_file)
    N²_end = diff(les_end.b[end]) ./ diff(les_end.zc)
    zf_end = (les_end.zc[1:end-1] .+ les_end.zc[2:end]) ./ 2
    below = findall(z -> z < -5, zf_end)
    z_base = zf_end[below[argmax(N²_end[below])]]
    κc_les, κu_les = nearest_value(les.zf, les.κc, z_base), nearest_value(les.zf, les.κu, z_base)
    push!(rows, @sprintf("| %s | LES | %.1f | %.2f | %.2e | %.2e | %.2f | | |", case_name, -z_base,
                         nearest_value(les.zf, les.Ri, z_base), κc_les, κu_les, κu_les / κc_les))

    fig = Figure(size = (1600, 650))
    titles = ["Ri = N² / S²", "κc (m²/s), tracers", "κu (m²/s), momentum", "Pr = κu / κc"]
    panels = [Axis(fig[1, i]; title = titles[i], ylabel = i == 1 ? "z (m)" : "", xscale = i == 1 ? identity : log10)
            for i in 1:4]
    depth_shown = -0.6case.depth
    plot_profile!(ax, x, z; kw...) = (keep = isfinite.(x) .& (x .> 0) .& (z .> depth_shown);
                                      any(keep) && lines!(ax, x[keep], z[keep]; kw...))
    plot_profile!(panels[1], clamp.(les.Ri, 1e-3, 10), les.zf; color = :black, linewidth = 3, label = "LES")
    plot_profile!(panels[2], les.κc, les.zf; color = :black, linewidth = 3)
    plot_profile!(panels[3], les.κu, les.zf; color = :black, linewidth = 3)
    plot_profile!(panels[4], les.κu ./ les.κc, les.zf; color = :black, linewidth = 3)

    colours = Makie.wong_colors()
    for (c, (label, variant, spacing)) in enumerate(diagnosed_columns)
        if startswith(variant, "nemo")
            file = nemo_tke_column(case_name, case, nemo_variants[variant]; spacing, label = variant)
        else
            file = run_column(case_name, case, variant, column_closures[variant], spacing)
        end
        col = column_diffusivities(column_profiles(file))
        style = spacing == 1 ? :solid : :dash
        plot_profile!(panels[1], clamp.(col.Ri, 1e-3, 10), col.zf; color = colours[c], linestyle = style, label)
        plot_profile!(panels[2], col.κc, col.zf; color = colours[c], linestyle = style)
        plot_profile!(panels[3], col.κu, col.zf; color = colours[c], linestyle = style)
        plot_profile!(panels[4], col.κu ./ col.κc, col.zf; color = colours[c], linestyle = style)
        κc_col, κu_col = nearest_value(col.zf, col.κc, z_base), nearest_value(col.zf, col.κu, z_base)
        push!(rows, @sprintf("| %s | %s | | %.2f | %.2e | %.2e | %.2f | %.1f | %.1f |", case_name, label,
                             nearest_value(col.zf, col.Ri, z_base), κc_col, κu_col, κu_col / κc_col,
                             κc_col / κc_les, κu_col / κu_les))
    end
    for ax in panels
        hlines!(ax, [z_base]; color = :gray, linestyle = :dot)
        ylims!(ax, depth_shown, 0)
    end
    axislegend(panels[1]; position = :rb, labelsize = 9)
    Label(fig[0, :], "$case_name: diffusivities from the LES fluxes (black) and the columns, median over the second half " *
                     "(dotted line: base of the LES mixed layer)"; fontsize = 15)
    save(joinpath(folder, "diagnose_$(case_name).png"), fig)
    return nothing
end

function diagnose_report(case_names, folder)
    rows = String[]
    for name in case_names
        diagnose(name, folder, rows)
    end
    io = IOBuffer()
    println(io, "# Diffusivities at the base of the mixed layer: LES fluxes against the columns\n")
    println(io, "At the depth of the LES's maximum N² at the end of the run; medians over the second half. κ_LES from the ",
                "resolved fluxes. Ratios > 1: the column mixes more than the LES at the base.\n")
    println(io, "| Case | Profile | base (m) | Ri | κc (m²/s) | κu (m²/s) | Pr | κc / κc_LES | κu / κu_LES |")
    println(io, "|---|---|---|---|---|---|---|---|---|")
    foreach(r -> println(io, r), rows)
    text = String(take!(io))
    print(text)
    write(joinpath(folder, "diagnose_fluxes.md"), text)
    return text
end

diagnose_report(diagnosed_cases, output_folder())
@info "Wrote $(joinpath(output_folder(), "diagnose_fluxes.md")) and diagnose_<case>.png"
