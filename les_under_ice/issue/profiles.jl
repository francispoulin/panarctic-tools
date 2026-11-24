# Profiles for the Oceananigans issue (Greg Wagner asked for profiles, 23 Nov): one clean figure, three cases (rows) ×
# buoyancy, u, v (columns) at the end of the run: the LES against default CATKE, adapted CATKE (rotation limit 0.1,
# Cˢ 0.991, half the surface TKE flux), k-ε and NEMO's TKE (wave terms off), all single columns at 4 m.
# CATKE columns from ice_variants.jl (Δt = 10 min), k-ε and NEMO from columns.jl / nemo_model.jl (Δt = 1 min); the time
# step changes CATKE's depth by a few % at most (ice_variants_dt1 vs dt10). Reads files only:
#
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/issue/profiles.jl")
#
# Writes output/profiles.png (and .pdf).

using Oceananigans, Oceananigans.Units
using CairoMakie, JLD2, Printf, Statistics

include(joinpath(@__DIR__, "..", "cases.jl"))

profile_cases = [("rotating_stress",      "Stress 0.41 N/m², rotating (75° N), 1 day"),
                 ("arctic_ice_drag",      "Ice drifting at 0.2 m/s, 4 days"),
                 ("arctic_ice_drag_melt", "Ice drifting at 0.2 m/s with melt, 4 days")]

profile_columns(name, folder) = [
    ("Default CATKE",  joinpath(folder, "ice_variants", "$(name)_stock_dz4_dt10.jld2"),                     "#d62728", :solid),
    ("Adapted CATKE",  joinpath(folder, "ice_variants", "$(name)_cap0.1_kappa041_flux_half_dz4_dt10.jld2"), "#1f77b4", :solid),
    ("k-ε",            joinpath(folder, "$(name)_column_k_epsilon_dz4.jld2"),                                "#9467bd", :dash),
    ("NEMO TKE",       joinpath(folder, "$(name)_column_nemo_no_waves_dz4.jld2"),                           "#2ca02c", :dashdot)]

function profiles_figure(folder)
    fig = Figure(size = (1350, 420 * length(profile_cases) + 80), fontsize = 15)
    Label(fig[0, 1:3], "Mean profiles at the end of the run: LES (black) against single columns at 4 m";
          fontsize = 19, font = :bold)
    for (row, (name, title)) in enumerate(profile_cases)
        case = case_parameters(name)
        les = les_profiles(joinpath(folder, name * "_les.jld2"))
        n_les = length(les.t)
        h_les = case_depth(name, case, les.b[n_les], les.zc)
        z_bottom = -min(case.depth, max(3h_les, 40.0))

        ax_b = Axis(fig[row, 1]; xlabel = row == length(profile_cases) ? "buoyancy (m s⁻²)" : "", ylabel = "z (m)")
        ax_u = Axis(fig[row, 2]; xlabel = row == length(profile_cases) ? "u (m s⁻¹)" : "")
        ax_v = Axis(fig[row, 3]; xlabel = row == length(profile_cases) ? "v (m s⁻¹)" : "")
        Label(fig[row, 1:3, Top()], title; fontsize = 16, padding = (0, 0, 8, 0), halign = :left, font = :bold)
        for ax in (ax_b, ax_u, ax_v)
            ylims!(ax, z_bottom, 0)
            hlines!(ax, [-h_les]; color = (:black, 0.35), linestyle = :dot)
        end
        hideydecorations!(ax_u; grid = false); hideydecorations!(ax_v; grid = false)

        # Initial buoyancy (grey), then the columns, then the LES on top
        first_col = nothing
        for (label, file, color, style) in profile_columns(name, folder)
            if !isfile(file)
                alt = replace(file, "nemo_no_waves" => "nemo")
                isfile(alt) ? (file = alt) : (@warn "$name: no $label column ($file)"; continue)
            end
            col = column_profiles(file)
            first_col === nothing && (first_col = col)
            n = argmin(abs.(col.t .- les.t[n_les]))
            lines!(ax_b, col.b[n], col.zc; color, linestyle = style, linewidth = 2.2, label)
            lines!(ax_u, col.u[n], col.zc; color, linestyle = style, linewidth = 2.2, label)
            lines!(ax_v, col.v[n], col.zc; color, linestyle = style, linewidth = 2.2, label)
        end
        first_col !== nothing &&
            lines!(ax_b, first_col.b[1], first_col.zc; color = (:grey, 0.6), linewidth = 1.5, label = "initial")
        lines!(ax_b, les.b[n_les], les.zc; color = :black, linewidth = 3.5, label = "LES (1 m)")
        lines!(ax_u, les.u[n_les], les.zc; color = :black, linewidth = 3.5, label = "LES (1 m)")
        lines!(ax_v, les.v[n_les], les.zc; color = :black, linewidth = 3.5, label = "LES (1 m)")
        text!(ax_b, 0.02, 0.03; text = @sprintf("LES mixed layer: %.0f m (dotted)", h_les), space = :relative,
              fontsize = 13, color = (:black, 0.7))
        row == 1 && axislegend(ax_b; position = :rb, framevisible = false, labelsize = 13)
    end
    colgap!(fig.layout, 18)
    save(joinpath(folder, "profiles.png"), fig; px_per_unit = 2)
    save(joinpath(folder, "profiles.pdf"), fig)
    return fig
end

profiles_figure(output_folder())
@info "Wrote $(joinpath(output_folder(), "profiles.png")) and .pdf"
