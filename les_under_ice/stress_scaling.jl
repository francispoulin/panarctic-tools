# Rung 3 of the validation ladder (claude/catke-les-plan.md): stress without waves into uniform N² = 1e-4 s⁻².
#   kato_phillips    no rotation: the LES and every column against Price's (1979) law h = 1.05 u★ (t/N)^½, as
#                    R(t) = h / (u★ (t/N)^½), which should be steady near 1.05
#   rotating_stress  f = 1.4e-4 s⁻¹: deepening should stop after about half an inertial period (~6 h)
# h: base of the mixed layer, the shallowest depth below 5 m where N² > N²/2, interpolated between faces (cases.jl).
# Pass (as rung 1, revised): |Δh| ≤ max(10 %, one cell) against the LES at the end, and b error ≤ 20 %.
# Writes output/stress_scaling.md and .png. After run_suite.jl with "rung3":
#
#     empty!(ARGS)
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/stress_scaling.jl")

using Oceananigans, Oceananigans.Units
using Oceananigans.OutputReaders: FieldTimeSeries
using CairoMakie, JLD2, Printf, Statistics

include(joinpath(@__DIR__, "cases.jl"))

function stress_coarse_grain(values, z_fine, zᶠ)
    return map(1:length(zᶠ)-1) do k
        inside = (z_fine .> zᶠ[k]) .& (z_fine .< zᶠ[k+1])
        any(inside) ? mean(values[inside]) : values[argmin(abs.(z_fine .- (zᶠ[k] + zᶠ[k+1]) / 2))]
    end
end

# Columns drawn in the figure (all are in the table)
shown_columns(label) = any(startswith(label, s) for s in ("stock_", "k_epsilon_", "nemo_no_waves_", "nemo_no_breaking_")) &&
                       (endswith(label, "dz1") || endswith(label, "nz50"))
column_colour(label) = startswith(label, "stock") ? :blue : startswith(label, "k_epsilon") ? :green :
                       startswith(label, "nemo_no_breaking") ? :orange : :red

function stress_report(folder)
    rows = String[]
    fig = Figure(size = (1500, 650))
    for (c, name) in enumerate(("kato_phillips", "rotating_stress"))
        case = case_parameters(name)
        N², u★ = stress_N², sqrt(abs(case.τx))
        les_file = joinpath(folder, name * "_les.jld2")
        isfile(les_file) || (@warn "No LES for $name: skipped"; continue)
        les = les_profiles(les_file)
        price_scale(t) = u★ * sqrt(t / sqrt(N²))
        ax = Axis(fig[1, c]; title = name == "kato_phillips" ? "Kato–Phillips: R = h / (u★ (t/N)^½)" :
                                     "Rotating stress: h (m)", xlabel = "t (h)")
        name == "kato_phillips" && hlines!(ax, [1.05]; color = :gray, linestyle = :dash)
        ns = findall(≥(1hour), les.t)
        h_les = [mixed_layer_base(les.b[n], les.zc, case.depth, N²; interpolate = true) for n in ns]
        y_les = name == "kato_phillips" ? h_les ./ price_scale.(les.t[ns]) : h_les
        lines!(ax, les.t[ns] ./ hour, y_les; color = :black, linewidth = 3, label = "LES")
        name == "kato_phillips" && push!(rows, @sprintf("| %s | LES | | %.1f | | R = %.2f (range after 6 h %.2f–%.2f) | |",
                                                        name, h_les[end], y_les[end],
                                                        extrema(y_les[les.t[ns] .≥ 6hours])...))
        files = sort(filter(f -> startswith(f, name * "_column_") && endswith(f, ".jld2"), readdir(folder)))
        for file in files
            label = replace(file, name * "_column_" => "", ".jld2" => "")
            col = column_profiles(joinpath(folder, file))
            nc = findall(≥(1hour), col.t)
            h_col = [mixed_layer_base(col.b[n], col.zc, case.depth, N²; interpolate = true) for n in nc]
            if shown_columns(label)
                y = name == "kato_phillips" ? h_col ./ price_scale.(col.t[nc]) : h_col
                lines!(ax, col.t[nc] ./ hour, y; color = column_colour(label),
                       linestyle = endswith(label, "nz50") ? :dot : :solid, label = label)
            end
            n_les = argmin(abs.(les.t .- col.t[end]))
            b_les = stress_coarse_grain(les.b[n_les], les.zc, col.zf)
            top = col.zc .> -0.75case.depth
            err = sqrt(sum(abs2, (col.b[end] .- b_les)[top]) / sum(abs2, (b_les .- col.b[1])[top]))
            h_end = mixed_layer_base(col.b[end], col.zc, case.depth, N²; interpolate = true)
            h_ref = mixed_layer_base(b_les, col.zc, case.depth, N²; interpolate = true)
            pass = abs(h_end - h_ref) ≤ max(0.1h_ref, cell_thickness_at(h_ref, col.zf)) && err ≤ 0.2
            R = name == "kato_phillips" ? @sprintf("R = %.2f", h_end / price_scale(col.t[end])) : ""
            push!(rows, @sprintf("| %s | %s | %.0f %% | %.1f | %.1f | %+.0f %% %s | %s |", name, label, 100err, h_end, h_ref,
                                 100(h_end - h_ref) / h_ref, R, pass ? "pass" : "**fail**"))
        end
        c == 1 && axislegend(ax; position = :rb, labelsize = 9)
    end
    Label(fig[0, :], "Rung 3: stress without waves (CATKE blue, k-ε green, NEMO base red, NEMO without breaking orange; dotted nz50)";
          fontsize = 16)
    save(joinpath(folder, "stress_scaling.png"), fig)

    io = IOBuffer()
    println(io, "# Rung 3: stress without waves into uniform N² = 1e-4 s⁻²\n")
    println(io, "u★ = 0.02 m/s. h: base of the mixed layer (N² > N²/2), interpolated between faces. Kato–Phillips: ",
                "Price (1979) R = h / (u★ (t/N)^½) ≈ 1.05. Pass: |Δh| ≤ max(10 %, one cell) and b error ≤ 20 % at 24 h.\n")
    println(io, "| Case | Column | b error | h (m) | h LES (m) | Δh | |")
    println(io, "|---|---|---|---|---|---|---|")
    foreach(r -> println(io, r), rows)
    text = String(take!(io))
    print(text)
    write(joinpath(folder, "stress_scaling.md"), text)
    return text
end

stress_report(output_folder())
@info "Wrote $(joinpath(output_folder(), "stress_scaling.md")) and stress_scaling.png"
