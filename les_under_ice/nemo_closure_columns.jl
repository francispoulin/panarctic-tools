# Check of NEMO's TKE scheme in Oceananigans (patches 084 / 085: NEMOTKEVerticalDiffusivity, a flavour of CATKE) against
# our standalone NEMO column (nemo_tke.jl, checked against zdftke.F90 and a numpy port) and the LES, before any 3-D run.
#
# Same cases, forcing and grids; both without NEMO's wave terms (nemo_no_waves), as under ice. The two differ only in
# discretisation: the Oceananigans closure keeps e at cell centres (NEMO at faces), imposes the surface value through the
# top half cell, and treats dissipation fully implicitly. Expected: |h_closure − h_port| ≤ max(10 %, one cell).
# Then Δt = 10 min (the 15 km model's) for the closure, with the reference parameters and with the "glorys12" preset
# (patch 087: NEMO 3.1's TKE defaults, surface TKE 3.75 u★², surface mixing length 0.4 m; nemo_tke_presets.jl).
#
# Columns: <case>_column_nemo_closure_<grid>.jld2 (Δt = 1 min), ..._nemo_closure_dt10_<grid>.jld2,
# ..._nemo_closure_glorys12_dt10_<grid>.jld2; the port's
# <case>_column_nemo_no_waves_<grid>.jld2 is run if missing. Existing files are reused: delete them to rerun.
# Grids: 4 m, nz50 (the pan-Arctic model's) and glorys50 (GLORYS12's levels): glorys12 on glorys50 is the column
# closest to what GLORYS itself does.
# Writes output/nemo_closure_columns.md and nemo_closure_columns.png. CPU, roughly an hour.
#
# Run in a new REPL in panarctic-061 after patches 085 and 087 (nemo_tke_vertical_diffusivity.jl, nemo_tke_presets.jl):
#     cd /mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-061 && julia --project=.
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/nemo_closure_columns.jl")

using Oceananigans, Oceananigans.Units
using CairoMakie, JLD2, Printf, Statistics

include(joinpath(@__DIR__, "cases.jl"))
include(joinpath(@__DIR__, "column_model.jl"))
include(joinpath(@__DIR__, "nemo_model.jl"))      # the standalone port: NEMOTKEParameters, nemo_tke_column, nemo_variants

let file = joinpath(dirname(Base.active_project()), "nemo_tke_presets.jl")
    isfile(file) || error("$file not found: apply patches 085 and 087 in this clone (or start Julia in panarctic-061)")
    include(file)    # nemo_tke_preset, nemo_tke_closure (which loads the closure into Oceananigans)
end

nemo_closure(preset = "reference") = Base.invokelatest(nemo_tke_closure, preset)

nc_cases = ["kato_phillips", "rotating_stress", "arctic_ice_drag", "arctic_ice_drag_brine", "arctic_ice_drag_melt",
            "arctic_brine_convection", convection_case_name(1e-5, 2.4e-7), "med_wind_med_cooling", "strong_wind"]
nc_spacings = (4, "nz50", "glorys50")   # glorys50: GLORYS12's own 50 levels (1 m top cell)

function nc_coarse_grain(values, z_fine, zᶠ)
    return map(1:length(zᶠ)-1) do k
        inside = (z_fine .> zᶠ[k]) .& (z_fine .< zᶠ[k+1])
        any(inside) ? mean(values[inside]) : values[argmin(abs.(z_fine .- (zᶠ[k] + zᶠ[k+1]) / 2))]
    end
end

nc_file(name, label, spacing) = joinpath(output_folder(), @sprintf("%s_column_%s_%s.jld2", name, label, spacing_label(spacing)))

function run_nemo_closure_columns(folder)
    rows = String[]
    passes, total = 0, 0
    for name in nc_cases, spacing in nc_spacings
        les_file = joinpath(folder, name * "_les.jld2")
        isfile(les_file) || (@warn "No LES for $name: skipped"; continue)
        case = case_parameters(name)

        port = nc_file(name, "nemo_no_waves", spacing)
        isfile(port) || nemo_tke_column(name, case, nemo_variants["nemo_no_waves"]; spacing, label = "nemo_no_waves")
        closure_1 = nc_file(name, "nemo_closure", spacing)
        isfile(closure_1) || run_column(name, case, "nemo_closure", nemo_closure(), spacing; Δt = 1minute, file = closure_1)
        closure_10 = nc_file(name, "nemo_closure_dt10", spacing)
        isfile(closure_10) || run_column(name, case, "nemo_closure_dt10", nemo_closure(), spacing; Δt = 10minutes,
                                         file = closure_10)
        glorys12 = nc_file(name, "nemo_closure_glorys12_dt10", spacing)
        isfile(glorys12) || run_column(name, case, "nemo_closure_glorys12_dt10", nemo_closure("glorys12"), spacing;
                                       Δt = 10minutes, file = glorys12)
        catke = joinpath(folder, "ice_variants", "$(name)_stock_$(spacing_label(spacing))_dt10.jld2")

        les = les_profiles(les_file)
        depth(file) = begin
            col = column_profiles(file)
            b_les = nc_coarse_grain(les.b[argmin(abs.(les.t .- col.t[end]))], les.zc, col.zf)
            (h_les = case_depth(name, case, b_les, col.zc), h = case_depth(name, case, col.b[end], col.zc),
             cell = col.zf[end] - col.zf[end-1])
        end
        p, c1, c10, g12 = depth(port), depth(closure_1), depth(closure_10), depth(glorys12)
        ck = isfile(catke) ? depth(catke) : nothing
        h_les = c1.h_les
        # The cell at the port's depth, for the tolerance
        zf = column_profiles(closure_1).zf
        k = clamp(searchsortedfirst(zf, -p.h), 2, length(zf))
        cell = zf[k] - zf[k-1]
        ok = abs(c1.h - p.h) ≤ max(0.1p.h, cell)
        passes += ok
        total += 1
        pct(h) = @sprintf("%+.0f", 100 * (h - h_les) / h_les)
        push!(rows, @sprintf("| %s | %s | %.1f | %s | %s | %s | %s | %s | %+.1f (cell %.1f) %s |", name,
                             spacing_label(spacing), h_les, pct(p.h), pct(c1.h), pct(c10.h), pct(g12.h),
                             ck === nothing ? "–" : pct(ck.h),
                             c1.h - p.h, cell, ok ? "✓" : "✗"))
    end

    io = IOBuffer()
    println(io, "# NEMO's TKE scheme in Oceananigans (patch 084/085) against the standalone NEMO column and the LES\n")
    println(io, "nemo_closure_columns.jl. Depth error against the LES (%). Port: nemo_tke.jl, no wave terms, Δt = 1 min. ",
                "Closure: NEMOTKEVerticalDiffusivity at Δt = 1 and 10 min; glorys12: the same with nemo_parameters = \"glorys12\" ",
                "(NEMO 3.1 defaults: surface TKE 3.75 u★², surface mixing length 0.4 m), Δt = 10 min. CATKE: stock, Δt = 10 min (ice_variants.jl). ",
                "Last column: closure − port (m) at Δt = 1 min, ✓ if within max(10 %, one cell).\n")
    println(io, "Closure within tolerance of the port: $passes of $total.\n")
    println(io, "| Case | grid | h LES (m) | port | closure Δt 1 min | closure Δt 10 min | glorys12 Δt 10 min | ",
                "CATKE Δt 10 min | closure − port |")
    println(io, "|---|---|---|---|---|---|---|---|---|")
    foreach(r -> println(io, r), rows)
    text = String(take!(io))
    print(text)
    write(joinpath(folder, "nemo_closure_columns.md"), text)

    # Profiles at the end: b, u and v for three cases at 4 m
    shown = ["rotating_stress", "arctic_ice_drag", "arctic_ice_drag_melt"]
    fig = Figure(size = (1350, 420 * length(shown) + 60), fontsize = 14)
    Label(fig[0, 1:3], "NEMO TKE: Oceananigans closure against the standalone port and the LES (4 m)"; fontsize = 17,
          font = :bold)
    for (row, name) in enumerate(shown)
        les_file = joinpath(folder, name * "_les.jld2")
        isfile(les_file) || continue
        les = les_profiles(les_file)
        n = length(les.t)
        axes = [Axis(fig[row, c]; xlabel = row == length(shown) ? label : "", ylabel = c == 1 ? "z (m)" : "",
                     title = c == 1 ? name : "") for (c, label) in enumerate(("b (m s⁻²)", "u (m s⁻¹)", "v (m s⁻¹)"))]
        h = case_depth(name, case_parameters(name), les.b[n], les.zc)
        for ax in axes
            ylims!(ax, -max(3h, 40), 0)
        end
        for (file, color, style, label) in ((nc_file(name, "nemo_no_waves", 4), :green, :solid, "port, Δt 1 min"),
                                            (nc_file(name, "nemo_closure", 4), :blue, :dash, "closure, Δt 1 min"),
                                            (nc_file(name, "nemo_closure_dt10", 4), :red, :dot, "closure, Δt 10 min"),
                                            (nc_file(name, "nemo_closure_glorys12_dt10", 4), :orange, :dashdot,
                                             "glorys12, Δt 10 min"))
            isfile(file) || continue
            col = column_profiles(file)
            m = argmin(abs.(col.t .- les.t[n]))
            for (ax, ϕ) in zip(axes, (col.b, col.u, col.v))
                lines!(ax, ϕ[m], col.zc; color, linestyle = style, linewidth = 2.2, label)
            end
        end
        for (ax, ϕ) in zip(axes, (les.b, les.u, les.v))
            lines!(ax, ϕ[n], les.zc; color = :black, linewidth = 3, label = "LES")
        end
        row == 1 && axislegend(axes[1]; position = :rb, framevisible = false, labelsize = 12)
    end
    save(joinpath(folder, "nemo_closure_columns.png"), fig; px_per_unit = 2)
    return text
end

run_nemo_closure_columns(output_folder())
@info "Wrote $(joinpath(output_folder(), "nemo_closure_columns.md")) and nemo_closure_columns.png"
