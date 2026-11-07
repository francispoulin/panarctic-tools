# Compare the open boundary conditions: net transport, unbalanced energy left in the basin, error against the
# reference solution, and the domain-mean free surface. Reads results/*.jls written by run_case.
#
#   include("analysis.jl")
#   compare(:adjustment)                 # writes figures/adjustment.png and figures/adjustment_maps.png
#   compare(:offset); compare(:steady); compare(:net_inflow)

using Serialization, Printf, CairoMakie
isdefined(Main, :PVInverter) || include("pv_inversion.jl")

load(name, experiment; shelf = false, dir = "results") =
    deserialize(joinpath(dir, "$(name)_$(experiment)$(shelf ? "_shelf" : "").jls"))

basin_inverter(P; shelf = false) = begin
    N = round(Int, P.L / P.Δ)
    x = ((1:N) .- 0.5) .* P.Δ
    d(xi) = shelf ? P.H - (P.H - P.shelf_depth) * (1 - tanh((xi - P.shelf_x) / P.shelf_w)) / 2 : P.H
    PVInverter(N, N, P.Δ, P.Δ, P.f, P.g, [d(x[i]) for i in 1:N, j in 1:N]; open_south = true)
end

# Reference fields restricted to the basin (the reference domain extends south of y = 0)
function restrict(s, P)
    N = round(Int, P.L / P.Δ)
    j₀ = round(Int, P.south_extension / P.Δ)
    return (u = s.u[:, j₀+1:j₀+N], v = s.v[:, j₀+1:j₀+N+1], η = s.η[:, j₀+1:j₀+N])
end

# Analytic steady vortex on the basin's C-grid (truth for :steady and :offset; for :net_inflow plus a uniform rise)
function analytic_state(P; rise = 0.0)
    N = round(Int, P.L / P.Δ)
    xc = ((1:N) .- 0.5) .* P.Δ; xf = (0:N) .* P.Δ
    ηv(x, y) = P.η₀ * exp(-((x - P.x₀)^2 + (y - P.y₀)^2) / P.R^2)
    u = [(P.g / P.f) * ηv(x, y) * 2 * (y - P.y₀) / P.R^2 for x in xf, y in xc]
    v = [-(P.g / P.f) * ηv(x, y) * 2 * (x - P.x₀) / P.R^2 for x in xc, y in xf]
    η = [ηv(x, y) + rise for x in xc, y in xc]
    return (; u, v, η)
end

function reference_states(experiment, P, times; shelf = false, dir = "results")
    if experiment === :adjustment
        ref = load("reference", experiment; shelf, dir)
        return [restrict(s, P) for s in ref.snapshots]
    elseif experiment === :net_inflow
        return [analytic_state(P; rise = P.Q₀ * P.T_Q / 2π * (1 - cos(2π * t / P.T_Q)) / P.L^2) for t in times]
    else   # :steady, :offset, :offset_ramp: the steady vortex
        s = analytic_state(P)
        return [s for _ in times]
    end
end

rms(a) = sqrt(sum(abs2, a) / length(a))

function diagnostics(r, refs, inv, P)
    n = length(r.snapshots)
    E′ = zeros(n); EB = zeros(n); errη = zeros(n); E′ref = zeros(n)
    for k in 1:n
        s = r.snapshots[k]
        v_ext = Float64.(s.v[:, 1])                     # the model's own normal velocity on the open faces
        d = split_energy(inv, Float64.(s.u), Float64.(s.v), Float64.(s.η); v_ext)
        E′[k], EB[k] = d.E′, d.EB
        ref = refs[min(k, length(refs))]
        errη[k] = rms(s.η .- ref.η)
        E′ref[k] = split_energy(inv, ref.u, ref.v, ref.η; v_ext = Float64.(ref.v[:, 1])).E′
    end
    return (; E′, EB, errη, E′ref)
end

const labels = Dict("flather" => "M1 Flather", "dirichlet" => "M2 Dirichlet", "davies" => "M3 Davies",
                    "flather_target" => "M4 Flather + target", "balanced_net" => "M5b balanced (net)",
                    "balanced_pv" => "M5a balanced (PV)", "balanced_filtered" => "M6 filtered C (2 d)",
                    "balanced_filtered_T12h" => "M6 filtered C (12 h)", "balanced_filtered_T5d" => "M6 filtered C (5 d)")
const colours = Dict("flather" => :steelblue, "dirichlet" => :firebrick, "davies" => :darkorange,
                     "flather_target" => :purple, "balanced_net" => :seagreen, "balanced_pv" => :black,
                     "balanced_filtered" => :deeppink, "balanced_filtered_T12h" => :hotpink,
                     "balanced_filtered_T5d" => :mediumvioletred)

const default_methods = ("flather", "dirichlet", "davies", "flather_target", "balanced_net", "balanced_pv", "balanced_filtered")

function compare(experiment; methods = default_methods, tag_suffix = "",
                 shelf = false, dir = "results", figdir = "figures")
    results = [load(m, experiment; shelf, dir) for m in methods if isfile(joinpath(dir, "$(m)_$(experiment)$(shelf ? "_shelf" : "").jls"))]
    isempty(results) && error("no results for $experiment")
    P = results[1].P
    inv = basin_inverter(P; shelf)
    refs = reference_states(experiment, P, results[1].times; shelf, dir)

    fig = Figure(size = (1400, 900), fontsize = 14)
    tag = string(experiment, shelf ? " (shelf)" : "")
    ax1 = Axis(fig[1, 1], title = "net transport into the basin, open side — $tag", xlabel = "day", ylabel = "Sv")
    ax2 = Axis(fig[1, 2], title = "unbalanced energy in the basin", xlabel = "day", ylabel = "J m³ kg⁻¹", yscale = log10)
    ax3 = Axis(fig[2, 1], title = "rms η error against the reference", xlabel = "day", ylabel = "m")
    ax4 = Axis(fig[2, 2], title = "domain-mean η", xlabel = "day", ylabel = "m")

    println("\n$tag: values at the last output")
    @printf("  %-22s %10s %12s %12s %10s\n", "method", "net [Sv]", "E′ [J m³/kg]", "rms Δη [cm]", "⟨η⟩ [cm]")
    for r in results
        t = r.times ./ 86400
        d = diagnostics(r, refs, inv, P)
        lab, col = labels[r.method], colours[r.method]
        lines!(ax1, t, r.net ./ 1e6; color = col, label = lab)
        lines!(ax2, t, max.(d.E′, 1e-6); color = col, label = lab)
        lines!(ax3, t, d.errη; color = col, label = lab)
        lines!(ax4, t, r.mean_η; color = col, label = lab)
        @printf("  %-22s %+10.3f %12.3e %12.2f %+10.2f\n", lab, r.net[end] / 1e6, d.E′[end], 100 * d.errη[end], 100 * r.mean_η[end])
        if r === results[1]
            lines!(ax1, t, r.net_ext ./ 1e6; color = :gray, linestyle = :dash, label = "external")
            if experiment === :adjustment
                lines!(ax2, t, max.(d.E′ref, 1e-6); color = :gray, linestyle = :dot, label = "reference")
                lines!(ax1, t[1:1], [NaN]; color = :gray, linestyle = :dot, label = "reference (E′ panel)")
            end
        end
    end
    Legend(fig[1:2, 3], ax1)
    mkpath(figdir)
    path = joinpath(figdir, "$(experiment)$(shelf ? "_shelf" : "")$(tag_suffix).png")
    save(path, fig)
    println("  saved $path")

    # Maps of η − η_ref at day 1 and at the end
    times_k = [findmin(abs.(results[1].times .- 86400))[2], length(results[1].times)]
    fig2 = Figure(size = (300 * length(results) + 100, 620), fontsize = 13)
    N = round(Int, P.L / P.Δ)
    x = ((1:N) .- 0.5) .* P.Δ ./ 1e3
    lim = 0.02
    for (row, k) in enumerate(times_k), (col, r) in enumerate(results)
        Δη = r.snapshots[k].η .- refs[min(k, length(refs))].η
        ax = Axis(fig2[row, col], title = @sprintf("%s, day %.1f", labels[r.method], r.times[k] / 86400),
                  aspect = DataAspect(), xlabel = row == 2 ? "x (km)" : "", ylabel = col == 1 ? "y (km)" : "")
        hm = heatmap!(ax, x, x, Δη; colormap = :balance, colorrange = (-lim, lim))
        col == length(results) && Colorbar(fig2[row, col+1], hm, label = "η − η_ref (m)")
    end
    path2 = joinpath(figdir, "$(experiment)$(shelf ? "_shelf" : "")$(tag_suffix)_maps.png")
    save(path2, fig2)
    println("  saved $path2")
    return nothing
end
