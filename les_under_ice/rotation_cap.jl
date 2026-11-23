# Step 2b of the CATKE recalibration plan (claude/catke-les-plan.md, 23 Nov): a rotation limit on CATKE's stable length.
#
# Diagnosis (diagnose_fluxes.jl, 23 Nov): with rotation, the LES momentum diffusivity falls from 1e-2 near the surface
# to 1e-4 at 30 m, while CATKE keeps ~5e-2 down to 25 m; under ice drag CATKE keeps κ ~ 1e-2 through the 40 m weakly
# stratified layer where the LES turbulence has died. Nothing in CATKE's length scale knows about f. In a weakly
# stratified, rotating boundary layer eddies are limited by rotation (Rossby scaling, h ~ u★/f), not by N.
#
# The cap: in stable_length_scaleᶜᶜᶠ / ᶜᶜᶜ (catke_mixing_length.jl l. 76, 92), ℓ = min(Cˢ d, √e/N) becomes
#     ℓ = min(Cˢ d, √e/N, Cᶠ √e/|f|)
# - the same form as the stratification length, so it is local (no u★, which under ice depends on the ice–ocean relative
#   velocity and is not available in these functions);
# - f = 0 or Cᶠ = 0 leaves CATKE exactly as it is (Kato–Phillips and every non-rotating case are untouched);
# - it limits the momentum, tracer and TKE lengths together (all three use the stable length).
# The cap is an @eval override inside Oceananigans (loaded once per process). Cᶠ and f are read from a Ref set before
# each run: fine for these CPU columns, not a GPU or 3-D implementation (that would add Cᶠ to CATKEMixingLength and
# read f from the model's Coriolis).
#
# Variants: stock (Cᶠ = 0), Cᶠ = 0.03, 0.1, 0.3, 1, and Cᶠ = 0.3 with prandtl10, on the targets and guards of
# variants.jl, 4 m and nz50 grids. Writes output/rotation_cap.md and output/rotation_cap.png; columns in
# output/rotation_cap/ (reused if present: delete the folder to rerun). CPU, under an hour.
#
#     empty!(ARGS)                       # "dz4" for the 4 m grid only
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/rotation_cap.jl")
#
# Run it in its own REPL (not one used for panarctic.jl): the override stays loaded. With Cᶠ = 0 it reproduces stock.

using Oceananigans, Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEEquation, CATKEMixingLength
using CairoMakie, JLD2, Printf, Statistics

include(joinpath(@__DIR__, "cases.jl"))
include(joinpath(@__DIR__, "column_model.jl"))

cap_spacings = (!isempty(ARGS) && ARGS[1] == "dz4") ? (4,) : (4, "nz50")
empty!(ARGS)

const TKEModule = Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities

# Once per process: redefining a method again invalidates compiled code (workflow.md)
if !isdefined(TKEModule, :rotation_cap)
    @eval TKEModule begin
        const rotation_cap = Ref((Cᶠ = 0.0, f = 0.0))

        @inline function rotation_mixing_length(w★, FT)
            Cᶠ, f = rotation_cap[]
            off = (Cᶠ == 0) | (f == 0)
            return ifelse(off, FT(Inf), FT(Cᶠ * w★ / abs(f)))
        end

        @inline function stable_length_scaleᶜᶜᶠ(i, j, k, grid, closure, e, velocities, tracers, buoyancy)
            Cˢ = closure.mixing_length.Cˢ
            Cᵇ = closure.mixing_length.Cᵇ
            d_up   = Cˢ * depthᶜᶜᶠ(i, j, k, grid)
            d_down = Cᵇ * height_above_bottomᶜᶜᶠ(i, j, k, grid)
            d = min(d_up, d_down)
            ℓᴺ = stratification_mixing_lengthᶜᶜᶠ(i, j, k, grid, closure, e, tracers, buoyancy)
            w★ = ℑzᵃᵃᶠ(i, j, k, grid, turbulent_velocityᶜᶜᶜ, closure, e)
            ℓᶠ = rotation_mixing_length(w★, eltype(grid))
            ℓ = min(d, ℓᴺ, ℓᶠ)
            ℓ = ifelse(isnan(ℓ), d, ℓ)
            return ℓ
        end

        @inline function stable_length_scaleᶜᶜᶜ(i, j, k, grid, closure, e, velocities, tracers, buoyancy)
            Cˢ = closure.mixing_length.Cˢ
            Cᵇ = closure.mixing_length.Cᵇ
            d_up   = Cˢ * depthᶜᶜᶜ(i, j, k, grid)
            d_down = Cᵇ * height_above_bottomᶜᶜᶜ(i, j, k, grid)
            d = min(d_up, d_down)
            ℓᴺ = stratification_mixing_lengthᶜᶜᶜ(i, j, k, grid, closure, e, tracers, buoyancy)
            w★ = turbulent_velocityᶜᶜᶜ(i, j, k, grid, closure, e)
            ℓᶠ = rotation_mixing_length(w★, eltype(grid))
            ℓ = min(d, ℓᴺ, ℓᶠ)
            ℓ = ifelse(isnan(ℓ), d, ℓ)
            return ℓ
        end
    end
end

cap_catke(; mixing_length = (;)) = CATKEVerticalDiffusivity(mixing_length = CATKEMixingLength(; mixing_length...),
                                                            turbulent_kinetic_energy_equation = CATKEEquation())
# name => (closure, Cᶠ)
cap_variants = [
    "stock"            => (cap_catke(), 0.0),
    "cap0.03"          => (cap_catke(), 0.03),
    "cap0.1"           => (cap_catke(), 0.1),
    "cap0.3"           => (cap_catke(), 0.3),
    "cap1"             => (cap_catke(), 1.0),
    "cap0.3_prandtl10" => (cap_catke(mixing_length = (; Cʰⁱc = 0.0242)), 0.3)]

cap_targets = ["rotating_stress", "arctic_ice_drag", "arctic_ice_drag_brine", "arctic_ice_drag_melt",
               "weak_wind_strong_cooling", "med_wind_med_cooling", "strong_wind_weak_cooling", "strong_wind"]
cap_guards  = ["kato_phillips", convection_case_name(1e-5, 2.4e-7), convection_case_name(4e-5, 4.8e-7),
               "arctic_brine_convection", "weak_wind_strong_cooling_stokes", "med_wind_med_cooling_stokes",
               "strong_wind_weak_cooling_stokes", "strong_wind_stokes"]

function cap_coarse_grain(values, z_fine, zᶠ)
    return map(1:length(zᶠ)-1) do k
        inside = (z_fine .> zᶠ[k]) .& (z_fine .< zᶠ[k+1])
        any(inside) ? mean(values[inside]) : values[argmin(abs.(z_fine .- (zᶠ[k] + zᶠ[k+1]) / 2))]
    end
end

function run_rotation_cap(folder)
    cap_folder = mkpath(joinpath(folder, "rotation_cap"))
    names = first.(cap_variants)
    cases_run = vcat(cap_targets, cap_guards)
    errors = Dict{Tuple{String, String, String}, Float64}()      # (case, variant, grid) => Δh / h_LES
    rows = String[]
    for name in cases_run, spacing in cap_spacings
        case = case_parameters(name)
        les_file = joinpath(folder, name * "_les.jld2")
        isfile(les_file) || (@warn "No LES for $name: skipped"; continue)
        les = les_profiles(les_file)
        grid = spacing_label(spacing)
        h_les = NaN
        cells = String[]
        for (variant, (closure, Cᶠ)) in cap_variants
            file = joinpath(cap_folder, "$(name)_$(variant)_$(grid).jld2")
            if !isfile(file)
                TKEModule.rotation_cap[] = (Cᶠ = Cᶠ, f = case.f)
                try
                    run_column(name, case, variant, closure, spacing; file)
                finally
                    TKEModule.rotation_cap[] = (Cᶠ = 0.0, f = 0.0)
                end
            end
            col = column_profiles(file)
            b_les = cap_coarse_grain(les.b[argmin(abs.(les.t .- col.t[end]))], les.zc, col.zf)
            h_les = case_depth(name, case, b_les, col.zc)
            h_col = case_depth(name, case, col.b[end], col.zc)
            errors[(name, variant, grid)] = (h_col - h_les) / h_les
            push!(cells, @sprintf("%+.0f", 100errors[(name, variant, grid)]))
        end
        role = name in cap_targets ? "target" : "guard"
        push!(rows, "| $name | $role | $grid | $(@sprintf("%.1f", h_les)) | f = $(case.f) | " * join(cells, " | ") * " |")
    end

    summary = map(names) do variant
        target = [abs(errors[k]) for k in keys(errors) if k[2] == variant && k[1] in cap_targets]
        guard  = [abs(errors[k]) for k in keys(errors) if k[2] == variant && k[1] in cap_guards]
        @sprintf("| %s | %.0f %% | %.0f %% | %.0f %% |", variant, 100mean(target), 100maximum(target), 100maximum(guard))
    end

    io = IOBuffer()
    println(io, "# CATKE with a rotation cap on the stable length: depth error against the LES (%)\n")
    println(io, "ℓ = min(Cˢ d, √e/N, Cᶠ √e/|f|) (rotation_cap.jl). Stock = Cᶠ 0. Targets and guards as in variants.jl.\n")
    println(io, "| Variant | mean abs, targets | worst, targets | worst, guards |")
    println(io, "|---|---|---|---|")
    foreach(r -> println(io, r), summary)
    println(io, "\n| Case | role | grid | h LES (m) | f | " * join(names, " | ") * " |")
    println(io, "|---|---|---|---|---|" * repeat("---|", length(names)))
    foreach(r -> println(io, r), rows)
    text = String(take!(io))
    print(text)
    write(joinpath(folder, "rotation_cap.md"), text)

    grid = spacing_label(first(cap_spacings))
    shown = [c for c in cases_run if haskey(errors, (c, "stock", grid))]
    values = [100errors[(c, v, grid)] for c in shown, v in names]
    fig = Figure(size = (1300, 60 + 42length(shown)))
    ax = Axis(fig[1, 1]; xticks = (1:length(names), names), yticks = (1:length(shown), shown), yreversed = true,
              xticklabelrotation = π / 8, title = "CATKE rotation cap, depth error (%) at $grid: red too deep, blue too shallow")
    heatmap!(ax, 1:length(names), 1:length(shown), permutedims(values); colormap = [:royalblue, :white, :firebrick],
             colorrange = (-40, 40))
    for i in eachindex(shown), j in eachindex(names)
        text!(ax, j, i; text = @sprintf("%+.0f", values[i, j]), align = (:center, :center), fontsize = 12)
    end
    hlines!(ax, [length(intersect(shown, cap_targets)) + 0.5]; color = :black, linewidth = 2)
    save(joinpath(folder, "rotation_cap.png"), fig)
    return text
end

run_rotation_cap(output_folder())
@info "Wrote $(joinpath(output_folder(), "rotation_cap.md")) and rotation_cap.png"
