# Step 2 of the CATKE recalibration plan (claude/catke-les-plan.md, 22 Nov): test the suspects directly.
#
# Hypothesis: under ice (rotating shear, no waves) CATKE entrains too much at the base of the mixed layer because, in
# stable stratification, (1) its tracer mixing stays strong (Prandtl number only 1 → 2.5, reached at Ri ≥ 1.27; NEMO's
# reaches 10), (2) its dissipation weakens (CˡᵒD 1.604 → CʰⁱD 0.579: a dissipation length ~3 times longer), and (3) it
# carries TKE down efficiently (Cˡᵒe = 7.863). Each variant changes one of these through CATKE's parameters only:
#   stock                  Oceananigans' defaults
#   prandtl10              Cʰⁱc = 0.0242: high-Ri Prandtl number Cʰⁱu / Cʰⁱc = 10, as NEMO's
#   prandtl4               Cʰⁱc = 0.0605: Prandtl number 4
#   flat_dissipation       CʰⁱD = CˡᵒD = 1.604: dissipation does not weaken with stratification
#   fast_transition        CRi⁰ = 0.1, CRiᵟ = 0.4: the high-Ri values reached by Ri = 0.5 instead of 1.27
#   weak_tke_transport     Cˡᵒe = 2.0 (from 7.863)
#   prandtl10_flat         prandtl10 and flat_dissipation
#   all                    prandtl10, flat_dissipation, fast_transition and weak_tke_transport
# Targets (no waves; CATKE too deep): rotating_stress, arctic_ice_drag[_brine, _melt], and the four no-wave wind cases.
# Guards (CATKE already right; a fix must not break them): kato_phillips, two convection cases, arctic_brine_convection,
# and the four wind cases with waves.
#
# Columns only, on the 4 m grid and the pan-Arctic nz50 grid (CPU; under an hour; files in output/variants/, reused if
# present: delete the folder to rerun). Writes output/variants.md (depth error per case and variant, summaries) and
# output/variants.png (the scorecard for the variants).
#
#     empty!(ARGS)                       # "dz4" for the 4 m grid only
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/variants.jl")

using Oceananigans, Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEEquation, CATKEMixingLength
using CairoMakie, JLD2, Printf, Statistics

include(joinpath(@__DIR__, "cases.jl"))
include(joinpath(@__DIR__, "column_model.jl"))

variant_spacings = (!isempty(ARGS) && ARGS[1] == "dz4") ? (4,) : (4, "nz50")
empty!(ARGS)

catke(; mixing_length = (;), equation = (;)) =
    CATKEVerticalDiffusivity(mixing_length = CATKEMixingLength(; mixing_length...),
                             turbulent_kinetic_energy_equation = CATKEEquation(; equation...))
prandtl10 = (; Cʰⁱc = 0.0242)
flat_D    = (; CʰⁱD = 1.604)
fast_Ri   = (; CRi⁰ = 0.1, CRiᵟ = 0.4)
weak_e    = (; Cˡᵒe = 2.0)
catke_variants = [
    "stock"              => catke(),
    "prandtl10"          => catke(mixing_length = prandtl10),
    "prandtl4"           => catke(mixing_length = (; Cʰⁱc = 0.0605)),
    "flat_dissipation"   => catke(equation = flat_D),
    "fast_transition"    => catke(mixing_length = fast_Ri),
    "weak_tke_transport" => catke(mixing_length = weak_e),
    "prandtl10_flat"     => catke(mixing_length = prandtl10, equation = flat_D),
    "all"                => catke(mixing_length = (; prandtl10..., fast_Ri..., weak_e...), equation = flat_D)]

target_cases = ["rotating_stress", "arctic_ice_drag", "arctic_ice_drag_brine", "arctic_ice_drag_melt",
                "weak_wind_strong_cooling", "med_wind_med_cooling", "strong_wind_weak_cooling", "strong_wind"]
guard_cases  = ["kato_phillips", convection_case_name(1e-5, 2.4e-7), convection_case_name(4e-5, 4.8e-7),
                "arctic_brine_convection", "weak_wind_strong_cooling_stokes", "med_wind_med_cooling_stokes",
                "strong_wind_weak_cooling_stokes", "strong_wind_stokes"]

# case_depth, threshold_depth and uniform_N² are in cases.jl (interpolated between levels since 23 Nov)

function variant_coarse_grain(values, z_fine, zᶠ)
    return map(1:length(zᶠ)-1) do k
        inside = (z_fine .> zᶠ[k]) .& (z_fine .< zᶠ[k+1])
        any(inside) ? mean(values[inside]) : values[argmin(abs.(z_fine .- (zᶠ[k] + zᶠ[k+1]) / 2))]
    end
end

function run_variants(folder)
    variant_folder = mkpath(joinpath(folder, "variants"))
    names = first.(catke_variants)
    cases_run = vcat(target_cases, guard_cases)
    errors = Dict{Tuple{String, String, String}, Float64}()      # (case, variant, grid) => Δh / h_LES
    rows = String[]
    for name in cases_run, spacing in variant_spacings
        case = case_parameters(name)
        les_file = joinpath(folder, name * "_les.jld2")
        isfile(les_file) || (@warn "No LES for $name: skipped"; continue)
        les = les_profiles(les_file)
        grid = spacing_label(spacing)
        h_les = NaN
        cells = String[]
        for (variant, closure) in catke_variants
            file = joinpath(variant_folder, "$(name)_$(variant)_$(grid).jld2")
            isfile(file) || run_column(name, case, variant, closure, spacing; file)
            col = column_profiles(file)
            b_les = variant_coarse_grain(les.b[argmin(abs.(les.t .- col.t[end]))], les.zc, col.zf)
            h_les = case_depth(name, case, b_les, col.zc)
            h_col = case_depth(name, case, col.b[end], col.zc)
            errors[(name, variant, grid)] = (h_col - h_les) / h_les
            push!(cells, @sprintf("%+.0f", 100errors[(name, variant, grid)]))
        end
        role = name in target_cases ? "target" : "guard"
        push!(rows, "| $name | $role | $grid | $(@sprintf("%.1f", h_les)) | " * join(cells, " | ") * " |")
    end

    # Summaries per variant: mean |error| over targets, worst |error| over guards, both grids together
    summary = map(names) do variant
        target = [abs(errors[k]) for k in keys(errors) if k[2] == variant && k[1] in target_cases]
        guard  = [abs(errors[k]) for k in keys(errors) if k[2] == variant && k[1] in guard_cases]
        @sprintf("| %s | %.0f %% | %.0f %% | %.0f %% |", variant, 100mean(target), 100maximum(target), 100maximum(guard))
    end

    io = IOBuffer()
    println(io, "# CATKE variants: mixed-layer depth error against the LES (%)\n")
    println(io, "Variants change one suspect at a time (variants.jl). Targets: no waves, CATKE too deep. Guards: cases ",
                "CATKE already gets right. Depth as in stress_scaling.jl (uniform N²) or compare.jl (others).\n")
    println(io, "| Variant | mean abs, targets | worst, targets | worst, guards |")
    println(io, "|---|---|---|---|")
    foreach(r -> println(io, r), summary)
    println(io, "\n| Case | role | grid | h LES (m) | " * join(names, " | ") * " |")
    println(io, "|---|---|---|---|" * repeat("---|", length(names)))
    foreach(r -> println(io, r), rows)
    text = String(take!(io))
    print(text)
    write(joinpath(folder, "variants.md"), text)

    # Figure: cases × variants on the 4 m grid, red too deep, blue too shallow
    grid = spacing_label(first(variant_spacings))
    shown = [c for c in cases_run if haskey(errors, (c, "stock", grid))]
    values = [100errors[(c, v, grid)] for c in shown, v in names]
    fig = Figure(size = (1500, 60 + 42length(shown)))
    ax = Axis(fig[1, 1]; xticks = (1:length(names), names), yticks = (1:length(shown), shown), yreversed = true,
              xticklabelrotation = π / 8, title = "CATKE variants, depth error (%) at $grid: red too deep, blue too shallow")
    heatmap!(ax, 1:length(names), 1:length(shown), permutedims(values); colormap = [:royalblue, :white, :firebrick],
             colorrange = (-40, 40))
    for i in eachindex(shown), j in eachindex(names)
        text!(ax, j, i; text = @sprintf("%+.0f", values[i, j]), align = (:center, :center), fontsize = 12)
    end
    hlines!(ax, [length(intersect(shown, target_cases)) + 0.5]; color = :black, linewidth = 2)
    save(joinpath(folder, "variants.png"), fig)
    return text
end

run_variants(output_folder())
@info "Wrote $(joinpath(output_folder(), "variants.md")) and variants.png"
