# The 21-case Arctic LES matrix (les_matrix.jl, 23-24 Nov) plus the earlier stress-driven cases, against five column
# schemes on three grids at the production time step. The baseline for the calibration, and the test of what CATKE's
# under-ice error follows: f (rotation), u★ (drift speed) or the surface buoyancy flux.
#
# Schemes:
#   catke_stock      Oceananigans' CATKE defaults
#   catke_adapted    best of ice_variants.jl (23 Nov): rotation cap Cᶠ = 0.1, Cˢ = 0.991 (κ = 0.41), Cᵂu★ halved
#   nemo_reference   NEMO's TKE closure (patches 085, 087, 088), NEMO 3.6-4.x reference namelist
#   nemo_glorys12    the same with GLORYS12's NEMO 3.1 defaults (surface TKE 3.75 u★², surface mixing length 0.4 m)
#   k_epsilon        Oceananigans' TKEDissipationVerticalDiffusivity, defaults
# Grids: 4 m, nz50 (the pan-Arctic model's), glorys50 (GLORYS12's levels). Δt = 10 min (second argument to change).
#
# Depth error against the LES at the end of the run and at half time: if the error grows between the two, the scheme
# drifts away as the layer keeps deepening (the LES is still adjusting at the end: h f/u★ still rising slowly).
# Files in output/matrix_columns/ (reused if present: delete to rerun); writes output/matrix_columns.md, the heatmap
# matrix_columns.png and the sweeps matrix_columns_sweeps.png. CPU, workers (Distributed, one thread each).
#
# In a new REPL in panarctic-061 after patches 085, 087, 088 (nemo_tke_presets.jl loads the NEMO closure):
#     empty!(ARGS); push!(ARGS, "24", "10")        # workers, Δt in minutes
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/matrix_columns.jl")

using Distributed

mc_workers    = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 24
mc_Δt_minutes = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 10
empty!(ARGS)

mc_folder = @__DIR__
mc_presets = joinpath(dirname(Base.active_project()), "nemo_tke_presets.jl")
isfile(mc_presets) || error("$mc_presets not found: apply patches 085, 087 and 088 in panarctic-061 and start Julia there")

existing_workers = workers() == [1] ? 0 : nworkers()
existing_workers < mc_workers &&
    addprocs(mc_workers - existing_workers; exeflags = "--project=$(dirname(Base.active_project()))",
             env = ["JULIA_NUM_THREADS" => "1"])

@everywhere begin
    using Logging
    using Oceananigans, Oceananigans.Units
    using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity, TKEDissipationVerticalDiffusivity
    using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEEquation, CATKEMixingLength
    include(joinpath($mc_folder, "cases.jl"))
    include(joinpath($mc_folder, "column_model.jl"))
    include(joinpath($mc_folder, "rotation_cap_override.jl"))
    include($mc_presets)
end

using CairoMakie, JLD2, Printf, Statistics

# Built on the worker from the scheme's name (plain data only through pmap)
@everywhere function mc_closure(scheme)
    scheme == "catke_stock"    && return CATKEVerticalDiffusivity(), 0.0
    scheme == "catke_adapted"  && return CATKEVerticalDiffusivity(mixing_length = CATKEMixingLength(Cˢ = 0.991),
                                              turbulent_kinetic_energy_equation = CATKEEquation(Cᵂu★ = 1.590)), 0.1
    scheme == "nemo_reference" && return nemo_tke_closure("reference"), 0.0
    scheme == "nemo_glorys12"  && return nemo_tke_closure("glorys12"), 0.0
    scheme == "k_epsilon"      && return TKEDissipationVerticalDiffusivity(), 0.0
    error("Unknown scheme $scheme")
end

@everywhere function mc_run(name, scheme, spacing, Δt_minutes, file)
    case = case_parameters(name)
    try
        closure, Cᶠ = mc_closure(scheme)
        set_rotation_cap!(Cᶠ, case.f)
        with_logger(NullLogger()) do
            run_column(name, case, scheme, closure, spacing; Δt = Δt_minutes * Oceananigans.Units.minute, file)
        end
        return file
    catch err
        return sprint(showerror, err)
    finally
        set_rotation_cap!(0, 0)
    end
end

mc_schemes = ["catke_stock", "catke_adapted", "nemo_reference", "nemo_glorys12", "k_epsilon"]
mc_spacings = (4, "nz50", "glorys50")
mc_earlier = ["kato_phillips", "rotating_stress", "arctic_ice_drag", "arctic_ice_drag_brine", "arctic_ice_drag_melt"]

mc_regime(name) = name == "kato_phillips" ? "no rotation (Kato–Phillips)" :
                  startswith(name, "rotating_stress") ? "constant stress" :
                  endswith(name, "_brine") ? "ice drag + brine" :
                  endswith(name, "_melt") ? "ice drag + melt" : "ice drag"

function mc_coarse_grain(values, z_fine, zᶠ)
    return map(1:length(zᶠ)-1) do k
        inside = (z_fine .> zᶠ[k]) .& (z_fine .< zᶠ[k+1])
        any(inside) ? mean(values[inside]) : values[argmin(abs.(z_fine .- (zᶠ[k] + zᶠ[k+1]) / 2))]
    end
end

# Depth error of a column against the LES at the column's time nearest t (relative, and the LES depth)
function mc_error(name, les, col, t)
    case = case_parameters(name)
    n = argmin(abs.(col.t .- t))
    b_les = mc_coarse_grain(les.b[argmin(abs.(les.t .- col.t[n]))], les.zc, col.zf)
    h_les = case_depth(name, case, b_les, col.zc)
    h_col = case_depth(name, case, col.b[n], col.zc)
    return (h_col - h_les) / h_les, h_les
end

function run_matrix_columns(folder, Δt_minutes)
    out = mkpath(joinpath(folder, "matrix_columns"))
    finished(name) = !(name in les_matrix) || isfile(joinpath(folder, "les_matrix", name * ".done"))
    names = [n for n in vcat(mc_earlier, les_matrix) if isfile(joinpath(folder, n * "_les.jld2")) && finished(n)]
    file_of(name, scheme, spacing) = joinpath(out, "$(name)_$(scheme)_$(spacing_label(spacing))_dt$(Δt_minutes).jld2")

    jobs = [(name, scheme, spacing, Δt_minutes, file_of(name, scheme, spacing))
            for name in names, spacing in mc_spacings, scheme in mc_schemes if !isfile(file_of(name, scheme, spacing))]

    # Fail fast: both NEMO presets and the adapted CATKE on a worker, one short case each (~1–2 min)
    for scheme in ("nemo_glorys12", "catke_adapted")
        check = remotecall_fetch(mc_run, first(workers()), "kato_phillips", scheme, 4, Δt_minutes,
                                 joinpath(out, "check_$(scheme)_dt$(Δt_minutes).jld2"))
        endswith(check, ".jld2") || error("$scheme failed on a worker:\n$check")
    end
    @info "Checks passed; $(length(names)) cases, $(length(jobs)) columns to run on $(nworkers()) workers (Δt = $Δt_minutes min)"

    results = pmap(job -> mc_run(job...), jobs)
    failures = [r for r in results if !endswith(r, ".jld2")]
    if !isempty(failures)
        @warn "$(length(failures)) columns failed; distinct errors:"
        foreach(m -> println("  ", first(m, 300)), unique(failures))
    end

    # Errors at the end and at half time
    final = Dict{Tuple{String, String, String}, Float64}()
    half = Dict{Tuple{String, String, String}, Float64}()
    h_les_end = Dict{Tuple{String, String}, Float64}()
    rows = String[]
    for name in names
        les = les_profiles(joinpath(folder, name * "_les.jld2"))    # read once per case (1 m LES files are large)
        for spacing in mc_spacings
        grid = spacing_label(spacing)
        cells = String[]
        for scheme in mc_schemes
            file = file_of(name, scheme, spacing)
            if !isfile(file)
                push!(cells, "–")
                continue
            end
            col = column_profiles(file)
            final[(name, scheme, grid)], h_les_end[(name, grid)] = mc_error(name, les, col, col.t[end])
            half[(name, scheme, grid)], _ = mc_error(name, les, col, col.t[end] / 2)
            push!(cells, @sprintf("%+.0f (%+.0f)", 100final[(name, scheme, grid)], 100half[(name, scheme, grid)]))
        end
        push!(rows, "| $name | $(mc_regime(name)) | $grid | $(@sprintf("%.1f", get(h_les_end, (name, grid), NaN))) | " *
                    join(cells, " | ") * " |")
        end
    end

    regimes = unique(mc_regime.(names))
    summary = String[]
    for spacing in mc_spacings, regime in regimes
        grid = spacing_label(spacing)
        cells = map(mc_schemes) do scheme
            e = [v for (k, v) in final if k[2] == scheme && k[3] == grid && mc_regime(k[1]) == regime]
            e½ = [v for (k, v) in half if k[2] == scheme && k[3] == grid && mc_regime(k[1]) == regime]
            isempty(e) ? "–" : @sprintf("%.0f (%+.0f; ½: %+.0f)", 100mean(abs.(e)), 100mean(e), 100mean(e½))
        end
        n = count(c -> mc_regime(c) == regime, names)
        push!(summary, "| $grid | $regime | $n | " * join(cells, " | ") * " |")
    end

    io = IOBuffer()
    println(io, "# LES matrix against five column schemes, three grids, Δt = $Δt_minutes min\n")
    println(io, "matrix_columns.jl. Mixed-layer depth error against the LES (%; > 0 too deep). catke_adapted: rotation cap ",
                "Cᶠ = 0.1, Cˢ = 0.991, Cᵂu★ halved; nemo_glorys12: NEMO 3.1 defaults (surface TKE 3.75 u★², surface ",
                "mixing length 0.4 m). Failed columns: $(length(failures)).\n")
    println(io, "## By regime\n")
    println(io, "Mean |error| at the end (mean signed error; ½: mean signed error at half time).\n")
    println(io, "| Grid | Regime | cases | " * join(mc_schemes, " | ") * " |")
    println(io, "|---|---|---|" * repeat("---|", length(mc_schemes)))
    foreach(r -> println(io, r), summary)
    println(io, "\n## By case\n")
    println(io, "Error at the end (at half time).\n")
    println(io, "| Case | regime | grid | h LES (m) | " * join(mc_schemes, " | ") * " |")
    println(io, "|---|---|---|---|" * repeat("---|", length(mc_schemes)))
    foreach(r -> println(io, r), rows)
    text = String(take!(io))
    print(text)
    write(joinpath(folder, "matrix_columns.md"), text)

    # Heatmaps: cases × schemes, one per grid
    fig = Figure(size = (2100, 1300))
    for (p, spacing) in enumerate(mc_spacings)
        grid = spacing_label(spacing)
        shown = [c for c in names if haskey(final, (c, "catke_stock", grid))]
        isempty(shown) && continue
        values = [get(final, (c, s, grid), NaN) * 100 for c in shown, s in mc_schemes]
        ax = Axis(fig[1, p]; xticks = (1:length(mc_schemes), mc_schemes), yticks = (1:length(shown), shown),
                  yreversed = true, xticklabelrotation = π / 5, title = "$grid, Δt = $Δt_minutes min: depth error (%)")
        heatmap!(ax, 1:length(mc_schemes), 1:length(shown), permutedims(values);
                 colormap = [:royalblue, :white, :firebrick], colorrange = (-60, 60), nan_color = :grey85)
        for i in eachindex(shown), j in eachindex(mc_schemes)
            isfinite(values[i, j]) && text!(ax, j, i; text = @sprintf("%+.0f", values[i, j]),
                                           align = (:center, :center), fontsize = 10)
        end
    end
    save(joinpath(folder, "matrix_columns.png"), fig)

    # Sweeps: error against f (drift 0.2 m/s) and against the ice drift speed (f = 1.4e-4), nz50 and glorys50
    f_sweep(flux) = [(flux == "" ? [("arctic_ice_drag_f0", 0.0)] : [])...,      # no f = 0 with brine or melt
                     ("arctic_ice_drag_f0.7e-4" * flux, 0.7e-4), ("arctic_ice_drag_f1.0e-4" * flux, 1.0e-4),
                     ("arctic_ice_drag" * flux, 1.4e-4)]
    u_sweep(flux) = [("arctic_ice_drag_u0.1" * flux, 0.1), ("arctic_ice_drag" * flux, 0.2),
                     ("arctic_ice_drag_u0.3" * flux, 0.3), ("arctic_ice_drag_u0.4" * flux, 0.4)]
    colours = Dict(zip(mc_schemes, [:blue, :deepskyblue, :red, :orange, :green]))
    styles = Dict("" => :solid, "_brine" => :dash, "_melt" => :dot)
    fig2 = Figure(size = (1500, 950))
    for (col, spacing) in enumerate(("nz50", "glorys50")), (row, (sweep, xlabel)) in
            enumerate(((f_sweep, "f (10⁻⁴ s⁻¹), ice drift 0.2 m/s"), (u_sweep, "ice drift (m/s), f = 1.4 × 10⁻⁴ s⁻¹")))
        ax = Axis(fig2[row, col]; xlabel, ylabel = "depth error (%)",
                  title = "$spacing: solid ice drag, dashed + brine, dotted + melt")
        hlines!(ax, [0]; color = :grey, linewidth = 1)
        for scheme in mc_schemes, flux in ("", "_brine", "_melt")
            points = [(x, 100final[(c, scheme, spacing)]) for (c, x) in sweep(flux) if haskey(final, (c, scheme, spacing))]
            length(points) < 2 && continue
            x = first.(points) .* (row == 1 ? 1e4 : 1)
            scatterlines!(ax, x, last.(points); color = colours[scheme], linestyle = styles[flux],
                          label = flux == "" ? scheme : nothing)
        end
        row == 1 && col == 1 && axislegend(ax; position = :lt, framevisible = false, labelsize = 11)
    end
    save(joinpath(folder, "matrix_columns_sweeps.png"), fig2)
    return text
end

run_matrix_columns(output_folder(), mc_Δt_minutes)
@info "Wrote $(joinpath(output_folder(), "matrix_columns.md")), matrix_columns.png and matrix_columns_sweeps.png"
