# Compare the single columns (columns.jl, nemo_tke.jl) with the LES (les.jl) of a case: buoyancy and velocity profiles
# at a quarter, half and the end of the run, and a table of errors. The LES is averaged onto each column's cells first (coarse-graining). Writes
# output/<case>_comparison.md and output/<case>_comparison.png. On a CPU, after both have run:
#
#     empty!(ARGS); push!(ARGS, "free_convection")            # a second argument "small" or "fine" as given to les.jl
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/compare.jl")
#
# Errors, for each column against the LES, at the end of the run:
#   b error     ‖b_column − b_LES‖ / ‖b_LES − b_initial‖ over the top 160 m, or the top 3/4 of a shallower domain (the
#               fraction of the change that is wrong). For seawater cases b is computed from T and S (TEOS-10).
#   h           boundary-layer depth: the first depth below 10 m where b is below its 10 m value by more than 1e-5 m s⁻²
#               (about the 0.03 kg/m³ density step of the mixed-layer depth in the pan-Arctic analysis)
# Step 0 passes if the stock CATKE columns' b errors are around 10 % and their depths within a few metres of the LES's,
# at 1 to 8 m, as Wagner et al. report for the calibration cases.

using Oceananigans, Oceananigans.Units
using Oceananigans.OutputReaders: FieldTimeSeries
using CairoMakie, JLD2, Printf, Statistics

include(joinpath(@__DIR__, "cases.jl"))

case_name = isempty(ARGS) ? "free_convection" : ARGS[1]
size_name = length(ARGS) ≥ 2 ? ARGS[2] : "default"
folder    = output_folder()

les_file = joinpath(folder, case_name * (size_name == "default" ? "" : "_" * size_name) * "_les.jld2")
isfile(les_file) || error("No LES output $les_file: run les.jl first")
les   = les_profiles(les_file)
z_les = les.zc
Δ_les = z_les[2] - z_les[1]
case  = case_parameters(case_name)
z_top_error = -min(160.0, 0.75 * case.depth)

time_index(times, t) = argmin(abs.(times .- t))

# The LES averaged over the cells of a column grid with faces zᶠ (nearest LES value where a column cell is finer)
function coarse_grain(values, z_fine, zᶠ)
    return map(1:length(zᶠ)-1) do k
        inside = (z_fine .> zᶠ[k]) .& (z_fine .< zᶠ[k+1])
        any(inside) ? mean(values[inside]) : values[argmin(abs.(z_fine .- (zᶠ[k] + zᶠ[k+1]) / 2))]
    end
end

# Reference at 10 m (as the pan-Arctic MLD criterion), not the top cell: the LES's top cell holds an unresolved surface
# layer (much denser than below in free convection), which made the depth jump with the coarse-graining
function boundary_layer_depth(b, z; Δb = 1e-5, z_ref = -10.0)
    reference = b[argmin(abs.(z .- z_ref))]
    k = findlast(k -> z[k] < z_ref && b[k] < reference - Δb, eachindex(b))
    return k === nothing ? -z[1] : -z[k]
end

column_files = sort(filter(f -> startswith(f, case_name * "_column_") && endswith(f, ".jld2"), readdir(folder)))
isempty(column_files) && error("No column output for $case_name in $folder: run columns.jl first")

t_stop = les.t[end]
times = [t_stop / 4, t_stop / 2, t_stop]
fig = Figure(size = (1500, 900))
axes_b = [Axis(fig[1, i]; title = "b at $(round(Int, t / hour)) h", xlabel = "b (m s⁻²)", ylabel = i == 1 ? "z (m)" : "")
          for (i, t) in enumerate(times)]
axes_u = [Axis(fig[2, i]; title = "u at $(round(Int, t / hour)) h", xlabel = "u (m/s)", ylabel = i == 1 ? "z (m)" : "")
          for (i, t) in enumerate(times)]
for (i, t) in enumerate(times)
    n = time_index(les.t, t)
    lines!(axes_b[i], les.b[n], z_les; color = :black, linewidth = 3, label = "LES ($(Δ_les) m)")
    lines!(axes_u[i], les.u[n], z_les; color = :black, linewidth = 3, label = "LES ($(Δ_les) m)")
end

rows = String[]
colors = Makie.wong_colors()
for (c, name) in enumerate(column_files)
    label = replace(name, case_name * "_column_" => "", ".jld2" => "")
    col = column_profiles(joinpath(folder, name))
    zc, zᶠ = col.zc, col.zf
    style = occursin("nemo", label) ? :dashdot : occursin("nz", label) ? :dot : occursin("no_surface_shear", label) ? :dash : :solid
    for (i, t) in enumerate(times)
        n = argmin(abs.(col.t .- t))
        lines!(axes_b[i], col.b[n], zc; color = colors[mod1(c, 7)], linestyle = style, label)
        lines!(axes_u[i], col.u[n], zc; color = colors[mod1(c, 7)], linestyle = style, label)
    end
    # Errors at the end of the run, over the top 160 m (or 3/4 of the domain)
    n_col, n_les = length(col.t), time_index(les.t, col.t[end])
    b_les = coarse_grain(les.b[n_les], z_les, zᶠ)
    b_col = col.b[n_col]
    b_ini = col.b[1]                  # the column's initial profile (saved at t = 0)
    top = zc .> z_top_error
    err = sqrt(sum(abs2, (b_col .- b_les)[top]) / sum(abs2, (b_les .- b_ini)[top]))
    h_les, h_col = boundary_layer_depth(b_les, zc), boundary_layer_depth(b_col, zc)
    push!(rows, @sprintf("| %s | %.0f %% | %.1f | %.1f | %+.1f |", label, 100err, h_col, h_les, h_col - h_les))
end
axislegend(axes_b[1]; position = :rb, labelsize = 9)
Label(fig[0, :], "$case_name: single columns against the LES"; fontsize = 18)
figure_file = joinpath(folder, case_name * "_comparison.png")
save(figure_file, fig)

io = IOBuffer()
println(io, "# $case_name: single columns against the LES\n")
println(io, "LES: $(les_file); at the end of the run ($(round(Int, les.t[end] / hour)) h). b error: ‖b_column − b_LES‖ / ",
            "‖b_LES − b_initial‖ over the top $(round(Int, -z_top_error)) m, the LES averaged onto the column's cells. h: first depth below 10 m where b is below its 10 m ",
            "value by more than 1e-5 m s⁻².\n")
println(io, "| Column | b error | h column (m) | h LES (m) | Δh (m) |")
println(io, "|---|---|---|---|---|")
foreach(row -> println(io, row), rows)
text = String(take!(io))
print(text)
write(joinpath(folder, case_name * "_comparison.md"), text)
@info "Wrote $(joinpath(folder, case_name * "_comparison.md")) and $figure_file"
