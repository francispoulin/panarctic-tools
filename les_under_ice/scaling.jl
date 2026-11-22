# Rung 1 of the validation ladder (claude/catke-les-plan.md): free convection into uniform stratification. For the five
# convection cases of cases.jl (convection_suite), checks
#   1. the LES against the entrainment law h² N² = 2(1 + 2A) Jb t: C = h² N² / (Jb t) should be about constant in time
#      and the same for every case (2 without entrainment; 2.8 for A = 0.2), and A from the LES's own buoyancy flux,
#      A = −min_z ⟨w b⟩ / Jb (resolved flux only: implicit LES, no subfilter flux, so A is slightly low);
#   2. every column (CATKE variants, k-ε, NEMO TKE; 1, 4, 8 m and the pan-Arctic grids) against the LES: h / h_LES in
#      time, and the b error and Δh at 24 h (as compare.jl).
# h here is the base of the mixed layer (where the law's h is): the shallowest depth below 5 m where N² exceeds half
# the background N², from 3 h on. Pairs of cases share Jb / N² but differ in N: dimensionally h = (Jb t / N²)^½ F(N t),
# so the pairs test whether F matters. Expected (numpy check of nemo_tke.jl, 21 Nov): NEMO's C ≈ 2.0 at 1 m (no
# entrainment: it convects only through enhanced diffusion), 2.2–2.5 on coarser grids. Writes output/scaling.md and output/scaling.png. After run_suite.jl with "convection":
#
#     empty!(ARGS)
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/scaling.jl")
#
# Pass criteria (set before looking, 21 Nov): LES C constant to ±10 % after 6 h and the same across cases to ±10 %;
# a scheme passes rung 1 if, in every case and on every grid, |Δh| ≤ 10 % of h_LES and the b error ≤ 20 % at 24 h.

using Oceananigans, Oceananigans.Units
using Oceananigans.OutputReaders: FieldTimeSeries
using CairoMakie, JLD2, Printf, Statistics

include(joinpath(@__DIR__, "cases.jl"))

folder = output_folder()

# The shallowest face below 5 m where N² exceeds half the background N² (searching downward). Not the maximum of N²:
# a scheme without entrainment (NEMO: a step from 0 to N²) has no maximum (21 Nov, numpy check of nemo_tke.jl)
function mixed_layer_base(b, zc, depth, N²_background)
    for k in length(b)-1:-1:1          # zc ascending: from the top down
        zf = (zc[k] + zc[k+1]) / 2
        -zf < 5 && continue
        -zf > 0.8depth && break
        (b[k+1] - b[k]) / (zc[k+1] - zc[k]) > 0.5N²_background && return -zf
    end
    return NaN
end

function boundary_layer_depth(b, z; Δb = 1e-5, z_ref = -10.0)
    reference = b[argmin(abs.(z .- z_ref))]
    k = findlast(k -> z[k] < z_ref && b[k] < reference - Δb, eachindex(b))
    return k === nothing ? -z[1] : -z[k]
end

function coarse_grain(values, z_fine, zᶠ)
    return map(1:length(zᶠ)-1) do k
        inside = (z_fine .> zᶠ[k]) .& (z_fine .< zᶠ[k+1])
        any(inside) ? mean(values[inside]) : values[argmin(abs.(z_fine .- (zᶠ[k] + zᶠ[k+1]) / 2))]
    end
end

entrainment_times(t) = findall(≥(3hours), t)

# In a function: no soft-scope ambiguity with globals left in the REPL by other scripts
function scaling_report(folder)
    summary_rows, column_rows = String[], String[]
    fig = Figure(size = (1600, 1000))
    ax_C = Axis(fig[1, 1]; title = "LES: C = h² N² / (Jb t)", xlabel = "t (h)", ylabel = "C")
    hlines!(ax_C, [2.0, 2.8]; color = :gray, linestyle = :dash)
    colors = Makie.wong_colors()

    for (c, (N², Jb)) in enumerate(convection_suite)
        name = convection_case_name(N², Jb)
        case = case_parameters(name)
        les_file = joinpath(folder, name * "_les.jld2")
        isfile(les_file) || (@warn "No LES for $name: skipped"; continue)
        les = les_profiles(les_file)
        ns  = entrainment_times(les.t)
        h_les = [mixed_layer_base(les.b[n], les.zc, case.depth, N²) for n in ns]
        C = h_les .^ 2 .* N² ./ (Jb .* les.t[ns])
        lines!(ax_C, les.t[ns] ./ hour, C; color = colors[c], label = @sprintf("N² = %.1e, Jb = %.1e", N², Jb))

        WB = FieldTimeSeries(les_file, "wb")
        late = findall(≥(les.t[end] / 2), WB.times)
        A = mean(-minimum(Array(interior(WB[n], 1, 1, :))) / Jb for n in late)
        late_C = C[les.t[ns] .≥ 6hours]
        push!(summary_rows, @sprintf("| %s | %.1e | %.1e | %.1f | %.2f | %.2f–%.2f | %.2f | %.2f |", name, N², Jb, h_les[end],
                                     C[end], minimum(late_C), maximum(late_C), A, 2(1 + 2A)))

        ax = Axis(fig[1 + (c - 1) ÷ 3, 2 + (c - 1) % 3]; title = @sprintf("N² = %.1e, Jb = %.1e: h / h_LES", N², Jb),
                  xlabel = "t (h)")
        hlines!(ax, [0.9, 1.0, 1.1]; color = :gray, linestyle = :dash)
        files = sort(filter(f -> startswith(f, name * "_column_") && endswith(f, ".jld2"), readdir(folder)))
        for (j, file) in enumerate(files)
            label = replace(file, name * "_column_" => "", ".jld2" => "")
            col = column_profiles(joinpath(folder, file))
            h_col = [mixed_layer_base(col.b[n], col.zc, case.depth, N²) for n in eachindex(col.t)]
            ratio = [h_col[i] / h_les[argmin(abs.(les.t[ns] .- col.t[i]))] for i in eachindex(col.t) if col.t[i] ≥ 3hours]
            colour = startswith(label, "nemo") ? :red : startswith(label, "k_epsilon") ? :green : :blue
            lines!(ax, col.t[col.t .≥ 3hours] ./ hour, ratio; color = colour, linestyle = occursin("nz", label) ? :dot : :solid)

            # At the end: b error and Δh (compare.jl's definitions) against the LES on the column's cells
            n_les = argmin(abs.(les.t .- col.t[end]))
            b_les = coarse_grain(les.b[n_les], les.zc, col.zf)
            b_col, b_ini = col.b[end], col.b[1]
            top = col.zc .> -min(160.0, 0.75case.depth)
            err = sqrt(sum(abs2, (b_col .- b_les)[top]) / sum(abs2, (b_les .- b_ini)[top]))
            h_end, h_les_end = mixed_layer_base(b_col, col.zc, case.depth, N²), mixed_layer_base(b_les, col.zc, case.depth, N²)
            pass = abs(h_end - h_les_end) ≤ 0.1h_les_end && err ≤ 0.2
            push!(column_rows, @sprintf("| %s | %s | %.0f %% | %.1f | %.1f | %+.0f %% | %s |", name, label, 100err, h_end,
                                        h_les_end, 100(h_end - h_les_end) / h_les_end, pass ? "pass" : "**fail**"))
        end
    end
    axislegend(ax_C; position = :rb, labelsize = 9)
    Label(fig[0, :], "Rung 1: free convection into uniform N². Columns: CATKE blue, NEMO red, k-ε green; dotted: pan-Arctic grids"; fontsize = 18)
    save(joinpath(folder, "scaling.png"), fig)

    io = IOBuffer()
    println(io, "# Rung 1: free convection into uniform stratification\n")
    println(io, "h: base of the mixed layer (shallowest depth below 5 m where N² > N²_background / 2). C = h² N² / (Jb t): 2 without entrainment, ",
                "2(1 + 2A) with entrainment ratio A. A from the LES's resolved ⟨w b⟩ over the second half of the run.\n")
    println(io, "| Case | N² | Jb | h LES 24 h (m) | C at 24 h | C range after 6 h | A | 2(1 + 2A) |")
    println(io, "|---|---|---|---|---|---|---|---|")
    foreach(r -> println(io, r), summary_rows)
    println(io, "\nColumns at the end (pass: |Δh| ≤ 10 % and b error ≤ 20 %):\n")
    println(io, "| Case | Column | b error | h (m) | h LES (m) | Δh | |")
    println(io, "|---|---|---|---|---|---|---|")
    foreach(r -> println(io, r), column_rows)
    text = String(take!(io))
    print(text)
    write(joinpath(folder, "scaling.md"), text)
    return text
end

scaling_report(folder)
@info "Wrote $(joinpath(folder, "scaling.md")) and scaling.png"
