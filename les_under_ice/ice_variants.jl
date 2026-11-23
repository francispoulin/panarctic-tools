# Greg Wagner's suggestions on PanArcticShenanigans #48 (23 Nov) together with the rotation cap, on four vertical grids.
#
# Greg: (1) reduce the surface TKE flux and the surface mixing-length coefficient Cˢ; first guess Cˢ such that the von
# Kármán constant κ = Cˢ (Cˡᵒu³ / CˡᵒD)^¼ is 0.41 (stock 1.131 ↔ 0.468; 0.41 ↔ Cˢ = 0.991); a toggle between
# open-ocean and ice-covered values later. (2) More vertical resolution: GLORYS has ~1 m at the top.
# rotation_cap.jl (23 Nov): Cᶠ = 0.1 improved or left unchanged every case and guard (worst guard 19 → 12 %), not
# enough under ice drag (+45 %) or melt (+80 %); Cᶠ = 0.03 overshoots and breaks the wave guards; 0.3 and 1 never act.
#
# Variants (mixing length, TKE equation, Cᶠ):
#   stock                          Oceananigans' defaults
#   kappa041                       Cˢ = 0.991 (κ = 0.41)
#   flux_half                      Cᵂu★ = 1.590 (surface TKE flux from u★ halved)
#   flux_zero                      Cᵂu★ = 0, CᵂwΔ = 0 (no surface TKE flux)
#   kappa041_flux_half             both
#   cap0.1                         rotation cap Cᶠ = 0.1
#   cap0.1_kappa041                cap and Cˢ
#   cap0.1_kappa041_flux_half      all three
# Grids: 4 m, nz50, nz100 (the pan-Arctic model's) and glorys50 (GLORYS12's levels; cases.jl column_faces).
# Targets and guards as in variants.jl.
#
# Runs the columns on worker processes (Distributed, one thread each); files in output/ice_variants/ (reused if present:
# delete the folder to rerun); writes output/ice_variants_dt<Δt>.md and .png. ~30 min on 24 workers.
# Time step in minutes as the second argument: 1 (as every column so far) or 10 (the 15 km model's).
#
#     empty!(ARGS); push!(ARGS, "24", "1")        # workers, Δt in minutes
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/ice_variants.jl")

using Distributed

iv_workers    = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 24
iv_Δt_minutes = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 1
empty!(ARGS)

iv_folder = @__DIR__
existing_workers = workers() == [1] ? 0 : nworkers()
existing_workers < iv_workers &&
    addprocs(iv_workers - existing_workers; exeflags = "--project=$(dirname(Base.active_project()))",
             env = ["JULIA_NUM_THREADS" => "1"])

@everywhere begin
    using Logging
    using Oceananigans, Oceananigans.Units
    using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
    using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEEquation, CATKEMixingLength
    include(joinpath($iv_folder, "cases.jl"))
    include(joinpath($iv_folder, "column_model.jl"))
    include(joinpath($iv_folder, "rotation_cap_override.jl"))
end

using CairoMakie, JLD2, Printf, Statistics

# One column on a worker, from plain data only (screen.jl, 23 Nov: the function sent to pmap must not use anything
# defined only on this process). Returns the file, or the error text.
@everywhere function ice_variant_run(name, mixing_length, equation, Cᶠ, spacing, Δt_minutes, file)
    case = case_parameters(name)
    closure = CATKEVerticalDiffusivity(mixing_length = CATKEMixingLength(; mixing_length...),
                                       turbulent_kinetic_energy_equation = CATKEEquation(; equation...))
    set_rotation_cap!(Cᶠ, case.f)
    try
        with_logger(NullLogger()) do
            run_column(name, case, basename(file), closure, spacing; Δt = Δt_minutes * Oceananigans.Units.minute, file)
        end
        return file
    catch err
        return sprint(showerror, err)
    finally
        set_rotation_cap!(0, 0)
    end
end

κ041      = (; Cˢ = 0.991)
flux_half = (; Cᵂu★ = 1.590)
iv_variants = [
    "stock"                     => ((;),   (;),                          0.0),
    "kappa041"                  => (κ041,  (;),                          0.0),
    "flux_half"                 => ((;),   flux_half,                    0.0),
    "flux_zero"                 => ((;),   (; Cᵂu★ = 0.0, CᵂwΔ = 0.0),   0.0),
    "kappa041_flux_half"        => (κ041,  flux_half,                    0.0),
    "cap0.1"                    => ((;),   (;),                          0.1),
    "cap0.1_kappa041"           => (κ041,  (;),                          0.1),
    "cap0.1_kappa041_flux_half" => (κ041,  flux_half,                    0.1)]

iv_targets = ["rotating_stress", "arctic_ice_drag", "arctic_ice_drag_brine", "arctic_ice_drag_melt",
              "weak_wind_strong_cooling", "med_wind_med_cooling", "strong_wind_weak_cooling", "strong_wind"]
iv_guards  = ["kato_phillips", convection_case_name(1e-5, 2.4e-7), convection_case_name(4e-5, 4.8e-7),
              "arctic_brine_convection", "weak_wind_strong_cooling_stokes", "med_wind_med_cooling_stokes",
              "strong_wind_weak_cooling_stokes", "strong_wind_stokes"]
iv_spacings = (4, "nz50", "nz100", "glorys50")

function iv_coarse_grain(values, z_fine, zᶠ)
    return map(1:length(zᶠ)-1) do k
        inside = (z_fine .> zᶠ[k]) .& (z_fine .< zᶠ[k+1])
        any(inside) ? mean(values[inside]) : values[argmin(abs.(z_fine .- (zᶠ[k] + zᶠ[k+1]) / 2))]
    end
end

function run_ice_variants(folder, Δt_minutes)
    iv_folder_out = mkpath(joinpath(folder, "ice_variants"))
    names = first.(iv_variants)
    cases_run = [c for c in vcat(iv_targets, iv_guards) if isfile(joinpath(folder, c * "_les.jld2"))]
    length(cases_run) < length(iv_targets) + length(iv_guards) && @warn "Cases without LES skipped"
    file_of(name, variant, spacing) =
        joinpath(iv_folder_out, "$(name)_$(variant)_$(spacing_label(spacing))_dt$(Δt_minutes).jld2")

    jobs = [(name, ml, eq, Cᶠ, spacing, Δt_minutes, file_of(name, variant, spacing))
            for name in cases_run, spacing in iv_spacings, (variant, (ml, eq, Cᶠ)) in iv_variants
            if !isfile(file_of(name, variant, spacing))]

    # Fail fast: one column with the default parameters on a worker (~1 min)
    check = remotecall_fetch(ice_variant_run, first(workers()), "kato_phillips", (;), (;), 0.0, 4, Δt_minutes,
                             joinpath(iv_folder_out, "check_dt$(Δt_minutes).jld2"))
    endswith(check, ".jld2") || error("Default column failed on a worker:\n$check")
    @info "Check passed; $(length(jobs)) columns to run on $(nworkers()) workers (Δt = $Δt_minutes min)"

    results = pmap(job -> ice_variant_run(job...), jobs)
    failures = [r for r in results if !endswith(r, ".jld2")]
    if !isempty(failures)
        @warn "$(length(failures)) columns failed; distinct errors:"
        foreach(m -> println("  ", first(m, 300)), unique(failures))
    end

    errors = Dict{Tuple{String, String, String}, Float64}()      # (case, variant, grid) => Δh / h_LES
    rows = String[]
    for name in cases_run, spacing in iv_spacings
        case = case_parameters(name)
        les = les_profiles(joinpath(folder, name * "_les.jld2"))
        grid = spacing_label(spacing)
        h_les = NaN
        cells = String[]
        for variant in names
            file = file_of(name, variant, spacing)
            if !isfile(file)
                push!(cells, "–")
                continue
            end
            col = column_profiles(file)
            b_les = iv_coarse_grain(les.b[argmin(abs.(les.t .- col.t[end]))], les.zc, col.zf)
            h_les = case_depth(name, case, b_les, col.zc)
            h_col = case_depth(name, case, col.b[end], col.zc)
            errors[(name, variant, grid)] = (h_col - h_les) / h_les
            push!(cells, @sprintf("%+.0f", 100errors[(name, variant, grid)]))
        end
        role = name in iv_targets ? "target" : "guard"
        push!(rows, "| $name | $role | $grid | $(@sprintf("%.1f", h_les)) | " * join(cells, " | ") * " |")
    end

    # Summaries per variant and grid
    summary = String[]
    for spacing in iv_spacings, variant in names
        grid = spacing_label(spacing)
        target = [abs(v) for (k, v) in errors if k[2] == variant && k[3] == grid && k[1] in iv_targets]
        guard  = [abs(v) for (k, v) in errors if k[2] == variant && k[3] == grid && k[1] in iv_guards]
        isempty(target) && continue
        push!(summary, @sprintf("| %s | %s | %.0f %% | %.0f %% | %s |", grid, variant, 100mean(target),
                                100maximum(target), isempty(guard) ? "–" : @sprintf("%.0f %%", 100maximum(guard))))
    end

    io = IOBuffer()
    println(io, "# CATKE: Greg's suggestions (#48) and the rotation cap, four grids, Δt = $Δt_minutes min\n")
    println(io, "Depth error against the LES (%), ice_variants.jl. kappa041: Cˢ = 0.991; flux_half: Cᵂu★ = 1.59; ",
                "flux_zero: no surface TKE flux; cap0.1: ℓ ≤ 0.1 √e/|f|. Failed columns: $(length(failures)).\n")
    println(io, "| Grid | Variant | mean abs, targets | worst, targets | worst, guards |")
    println(io, "|---|---|---|---|---|")
    foreach(r -> println(io, r), summary)
    println(io, "\n| Case | role | grid | h LES (m) | " * join(names, " | ") * " |")
    println(io, "|---|---|---|---|" * repeat("---|", length(names)))
    foreach(r -> println(io, r), rows)
    text = String(take!(io))
    print(text)
    write(joinpath(folder, "ice_variants_dt$(Δt_minutes).md"), text)

    # Figure: one panel per grid, cases × variants
    fig = Figure(size = (2000, 1500))
    for (p, spacing) in enumerate(iv_spacings)
        grid = spacing_label(spacing)
        shown = [c for c in cases_run if haskey(errors, (c, "stock", grid))]
        isempty(shown) && continue
        values = [get(errors, (c, v, grid), NaN) * 100 for c in shown, v in names]
        ax = Axis(fig[(p - 1) ÷ 2 + 1, (p - 1) % 2 + 1]; xticks = (1:length(names), names),
                  yticks = (1:length(shown), shown), yreversed = true, xticklabelrotation = π / 6,
                  title = "$grid, Δt = $Δt_minutes min: depth error (%), red too deep")
        heatmap!(ax, 1:length(names), 1:length(shown), permutedims(values);
                 colormap = [:royalblue, :white, :firebrick], colorrange = (-40, 40), nan_color = :grey85)
        for i in eachindex(shown), j in eachindex(names)
            isfinite(values[i, j]) && text!(ax, j, i; text = @sprintf("%+.0f", values[i, j]),
                                           align = (:center, :center), fontsize = 10)
        end
        hlines!(ax, [length(intersect(shown, iv_targets)) + 0.5]; color = :black, linewidth = 2)
    end
    save(joinpath(folder, "ice_variants_dt$(Δt_minutes).png"), fig)
    return text
end

run_ice_variants(output_folder(), iv_Δt_minutes)
@info "Wrote $(joinpath(output_folder(), "ice_variants_dt$(iv_Δt_minutes).md")) and the .png"
