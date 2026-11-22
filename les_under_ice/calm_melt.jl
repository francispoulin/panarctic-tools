# Rung 2 of the validation ladder (claude/catke-les-plan.md): melt under calm conditions. Ice at rest over the Arctic
# profile, freshwater of 2 cm/day melt (Js = +5.2e-6 g/kg m/s), no stress, 4 days. No turbulence source, so no LES is
# needed: the freshwater should stay in the top cell (molecular diffusion spreads it ~2 cm in 4 days). The test is
# whether a scheme mixes when it should not.
#
# Runs the CATKE / k-ε columns and the NEMO columns for arctic_calm_melt (CPU, a few minutes; no GPU), then measures,
# from S − S_initial at the end:
#   outside top cell   fraction of the added freshwater below the top cell (truth ≈ 0)
#   d90                depth containing 90 % of the freshwater (truth: within the top cell)
#   budget             ∫ (S_initial − S) dz / (Js t): should be 1 (conservation check of the column itself)
# Writes output/calm_melt.md and output/calm_melt.png. In a REPL (CPU is enough):
#
#     empty!(ARGS)
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/calm_melt.jl")

using Oceananigans, Oceananigans.Units
using CairoMakie, JLD2, Printf

include(joinpath(@__DIR__, "cases.jl"))

calm_case_name = "arctic_calm_melt"
for script in ("columns.jl", "nemo_tke.jl")
    empty!(ARGS); push!(ARGS, calm_case_name)
    Base.include(Main, joinpath(@__DIR__, script))
end
empty!(ARGS)

function calm_melt_report(folder, name)
    case = case_parameters(name)
    rows = String[]
    fig = Figure(size = (1200, 700))
    ax = Axis(fig[1, 1]; title = "S − S_initial at $(round(Int, case.stop_time / day)) days", xlabel = "ΔS (g/kg)",
              ylabel = "z (m)")
    colours = Makie.wong_colors()
    files = sort(filter(f -> startswith(f, name * "_column_") && endswith(f, ".jld2"), readdir(folder)))
    for (c, file) in enumerate(files)
        label = replace(file, name * "_column_" => "", ".jld2" => "")
        col = column_profiles(joinpath(folder, file))
        Δz = diff(col.zf)                                    # cell thicknesses, bottom up
        fresh = (col.S[1] .- col.S[end]) .* Δz               # freshwater (as salt deficit) per cell
        total = sum(fresh)
        budget = total / (case.Js * col.t[end])
        top_cell = Δz[end]
        outside = sum(fresh[1:end-1]) / total
        # d90: from the surface down, the depth by which 90 % of the deficit is reached
        cumulative = cumsum(reverse(fresh)) ./ total          # top down
        k90 = findfirst(≥(0.9), cumulative)
        d90 = k90 === nothing ? NaN : -reverse(col.zf)[k90 + 1]
        push!(rows, @sprintf("| %s | %.1f | %.1f %% | %.1f | %.3f |", label, top_cell, 100outside, d90, budget))
        if endswith(label, "dz1") || endswith(label, "nz50")
            lines!(ax, col.S[end] .- col.S[1], col.zc; color = colours[mod1(c, 7)],
                   linestyle = endswith(label, "nz50") ? :dot : :solid, label = label)
        end
    end
    ylims!(ax, -40, 0)
    Legend(fig[1, 2], ax; labelsize = 9)
    save(joinpath(folder, "calm_melt.png"), fig)

    io = IOBuffer()
    println(io, "# Rung 2: melt under calm conditions ($name)\n")
    println(io, "Freshwater of 2 cm/day melt for 4 days, no stress. Truth: all of it in the top cell (molecular ",
                "diffusion ~2 cm). Outside top cell: fraction of the freshwater below the top cell. d90: depth holding ",
                "90 % of it. Budget: freshwater in the column / freshwater put in (1 = conserved).\n")
    println(io, "| Column | top cell (m) | outside top cell | d90 (m) | budget |")
    println(io, "|---|---|---|---|---|")
    foreach(r -> println(io, r), rows)
    text = String(take!(io))
    print(text)
    write(joinpath(folder, "calm_melt.md"), text)
    return text
end

calm_melt_report(output_folder(), calm_case_name)
@info "Wrote $(joinpath(output_folder(), "calm_melt.md")) and calm_melt.png"
