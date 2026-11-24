# Two diagnostics of rotation in the LES (23 Nov, after Greg Wagner's reply on the Oceananigans issue). LES files only.
#
# 1. Stress–shear angle. Every column scheme assumes a down-gradient (scalar) eddy viscosity: the turbulent stress
#    τ = −(⟨w u⟩, ⟨w v⟩) is parallel to the shear S = (∂z U, ∂z V). Greg's suggestion adds a cross-gradient term,
#        u′w′ = −K ∂z U + L² f ∂z V   (and v′w′ = −K ∂z V − L² f ∂z U),   i.e.   τ = K S + L² f ẑ × S,
#    so the stress is rotated from the shear by θ with tan θ = L² f / K. From the LES, on interior faces at every output
#    time of the second half of the run:
#        K∥ = τ · S / |S|²,   K⊥ = (ẑ × S) · τ / |S|² = (S_x τ_y − S_y τ_x) / |S|²,   θ = atan(K⊥, K∥),   L = √(|K⊥| / |f|)
#    θ > 0: stress rotated anticlockwise (to the left) of the shear. Values are kept where |S| and |τ| exceed 5 % of their
#    maximum at that time, and below 3 m (the resolved fluxes of an implicit LES are too small near the surface).
#    Binned by depth relative to the LES mixed-layer depth h(t); medians in the interior (0.3–0.8 h) and at the base
#    (0.8–1.1 h), as diagnose_fluxes.jl.
# 2. Depth against rotating boundary-layer scalings, over the second half of the run:
#        neutral rotating (Ekman) layer   h ≈ c u★/|f|, c ≈ 0.25–0.4 (Rossby and Montgomery; McPhee under ice)
#        stratified, rotating             h ≈ 1.7 u★/√(N f)  (Pollard, Rhines and Thompson 1973), N below the layer
#        stratified, no rotation          h ≈ 1.05 u★ (t/N)^½ (Price 1979)
#    u★ = √|τ| for a constant stress; under drifting ice u★ = √C_d |u_ice − u_water(top)| from the LES itself.
#
# Cases: every stress-driven LES in output/ (kato_phillips, rotating_stress*, arctic_ice_drag*), so it covers the
# matrix as it finishes. Writes output/rotation_diagnostics.md, stress_shear_angle.png and depth_scaling.png.
#
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/rotation_diagnostics.jl")

using Oceananigans, Oceananigans.Units
using Oceananigans.OutputReaders: FieldTimeSeries
using CairoMakie, JLD2, Printf, Statistics

include(joinpath(@__DIR__, "cases.jl"))

rd_nanmedian(x) = (y = filter(isfinite, x); isempty(y) ? NaN : median(y))
rd_quantile(x, q) = (y = filter(isfinite, x); isempty(y) ? NaN : quantile(y, q))
rd_series(file, name) = (fts = FieldTimeSeries(file, name); [vec(Array(interior(fts[n], 1, 1, :))) for n in eachindex(fts.times)])
function rd_interior_faces(values, Nz)
    length(values) == Nz     && return (values[1:end-1] .+ values[2:end]) ./ 2
    length(values) == Nz + 1 && return values[2:end-1]
    error("Flux profile of length $(length(values)) for Nz = $Nz")
end

# Stress-driven cases with an LES; a matrix case only once les_matrix.jl has marked it done (not while it is written)
rd_finished(name, folder) = !(name in les_matrix) || isfile(joinpath(folder, "les_matrix", name * ".done"))
rd_cases(folder) = sort([name for name in keys(cases)
                         if (name == "kato_phillips" || startswith(name, "rotating_stress") || startswith(name, "arctic_ice_drag")) &&
                            isfile(joinpath(folder, name * "_les.jld2")) && rd_finished(name, folder)])
rd_flux_label(name) = occursin("brine", name) ? "brine" : occursin("melt", name) ? "melt" : "stress only"

# Initial N² just below a depth h (mean over 10 m), from the case's initial profile
function rd_N²_below(case, h)
    z = range(-h - 10, -h; length = 21)
    b = case.seawater ? [seawater_buoyancy(case.T_init(zi), case.S_init(zi)) for zi in z] : case.b_init.(z)
    return max((b[end] - b[1]) / 10, 1e-8)
end

function rotation_diagnostics(name, folder)
    case = case_parameters(name)
    file = joinpath(folder, name * "_les.jld2")
    les = les_profiles(file)
    wu, wv = rd_series(file, "wu"), rd_series(file, "wv")
    Nz = length(les.zc)
    zf = (les.zc[1:end-1] .+ les.zc[2:end]) ./ 2
    Δz = diff(les.zc)
    late = findall(≥(les.t[end] / 2), les.t)
    f = case.f

    # 1. Angle, K∥, K⊥ binned by z / h
    bins = 0:0.1:1.2
    θ_bins = [Float64[] for _ in 1:length(bins)-1]
    θ_interior, θ_base, L_interior, ratio_interior = Float64[], Float64[], Float64[], Float64[]
    hs, u★s = Float64[], Float64[]
    for n in late
        h = case_depth(name, case, les.b[n], les.zc)
        push!(hs, h)
        Sx, Sy = diff(les.u[n]) ./ Δz, diff(les.v[n]) ./ Δz
        τx, τy = -rd_interior_faces(wu[n], Nz), -rd_interior_faces(wv[n], Nz)
        S, τ = hypot.(Sx, Sy), hypot.(τx, τy)
        keep = (S .> 0.05maximum(S)) .& (τ .> 0.05maximum(τ)) .& (zf .< -3)
        K∥ = (τx .* Sx .+ τy .* Sy) ./ S .^ 2
        K⊥ = (Sx .* τy .- Sy .* τx) ./ S .^ 2
        θ = atand.(K⊥, K∥)
        for k in findall(keep)
            ζ = -zf[k] / h
            b = searchsortedlast(bins, ζ)
            1 ≤ b < length(bins) && push!(θ_bins[b], θ[k])
            if 0.3 ≤ ζ ≤ 0.8
                push!(θ_interior, θ[k]); push!(ratio_interior, K⊥[k] / K∥[k])
                f != 0 && push!(L_interior, sqrt(abs(K⊥[k]) / abs(f)))
            elseif 0.8 < ζ ≤ 1.1
                push!(θ_base, θ[k])
            end
        end
        # u★: constant stress, or drag on the ice velocity relative to the water at the top
        push!(u★s, case.ice_drag ? sqrt(case.Cd) * hypot(case.uᵢ - les.u[n][end], case.vᵢ - les.v[n][end]) :
                                   sqrt(abs(case.τx)))
    end
    angle_profile = [(ζ = (bins[b] + bins[b+1]) / 2, median = rd_nanmedian(θ_bins[b]),
                      q25 = rd_quantile(θ_bins[b], 0.25), q75 = rd_quantile(θ_bins[b], 0.75), n = length(θ_bins[b]))
                     for b in 1:length(bins)-1]

    # 2. Depth against the scalings (medians over the second half)
    h, u★ = median(hs), median(u★s)
    t_mid = median(les.t[late])
    N = sqrt(rd_N²_below(case, h))
    ekman = f == 0 ? NaN : u★ / abs(f)
    prt   = f == 0 ? NaN : 1.7u★ / sqrt(N * abs(f))
    price = 1.05u★ * sqrt(t_mid / N)
    # h(t) over the whole run, for the figure
    h_t = [case_depth(name, case, les.b[n], les.zc) for n in eachindex(les.t)]
    u★_t = case.ice_drag ? [sqrt(case.Cd) * hypot(case.uᵢ - les.u[n][end], case.vᵢ - les.v[n][end]) for n in eachindex(les.t)] :
                           fill(sqrt(abs(case.τx)), length(les.t))

    return (; name, f, flux = rd_flux_label(name), angle_profile,
              θ_interior = rd_nanmedian(θ_interior), θ_interior_q = (rd_quantile(θ_interior, 0.25), rd_quantile(θ_interior, 0.75)),
              θ_base = rd_nanmedian(θ_base), ratio_interior = rd_nanmedian(ratio_interior),
              L_interior = rd_nanmedian(L_interior), h, u★, N, ekman, prt, price,
              t = les.t, h_t, u★_t)
end

function run_rotation_diagnostics(folder)
    names = rd_cases(folder)
    @info "Rotation diagnostics for $(length(names)) LES: $(join(names, ", "))"
    results = [rotation_diagnostics(name, folder) for name in names]

    io = IOBuffer()
    println(io, "# Rotation in the LES: stress–shear angle and depth scalings\n")
    println(io, "rotation_diagnostics.jl. Second half of each run. θ: angle of the turbulent stress from the shear ",
                "(0 = down-gradient; > 0 to the left); K⊥/K∥ = tan θ; L = √(|K⊥|/|f|), the length in Greg's L² f term. ",
                "Interior 0.3–0.8 h, base 0.8–1.1 h.\n")
    println(io, "## 1. Stress–shear angle\n")
    println(io, "| Case | f (s⁻¹) | surface flux | θ interior (25–75 %) | θ base | K⊥/K∥ interior | L interior (m) | h (m) |")
    println(io, "|---|---|---|---|---|---|---|---|")
    for r in results
        println(io, @sprintf("| %s | %.1e | %s | %+.0f° (%+.0f° to %+.0f°) | %+.0f° | %+.2f | %s | %.1f |", r.name, r.f, r.flux,
                             r.θ_interior, r.θ_interior_q..., r.θ_base, r.ratio_interior,
                             isfinite(r.L_interior) ? @sprintf("%.1f", r.L_interior) : "–", r.h))
    end
    println(io, "\n## 2. Depth against rotating boundary-layer scalings\n")
    println(io, "h, u★ (from the LES), N below the layer (initial profile): medians over the second half. ",
                "Expected h f/u★ ≈ 0.25–0.4 for a neutral rotating layer, h/h_PRT ≈ 1 for a stratified rotating one.\n")
    println(io, "| Case | f (s⁻¹) | surface flux | h (m) | u★ (m/s) | N² below (s⁻²) | h f/u★ | h / h_PRT | h / h_Price |")
    println(io, "|---|---|---|---|---|---|---|---|---|")
    for r in results
        println(io, @sprintf("| %s | %.1e | %s | %.1f | %.4f | %.1e | %s | %s | %.2f |", r.name, r.f, r.flux, r.h, r.u★, r.N^2,
                             isfinite(r.ekman) ? @sprintf("%.2f", r.h / r.ekman) : "–",
                             isfinite(r.prt) ? @sprintf("%.2f", r.h / r.prt) : "–", r.h / r.price))
    end
    text = String(take!(io))
    print(text)
    write(joinpath(folder, "rotation_diagnostics.md"), text)

    # Figure 1: θ against depth / h, one panel per surface flux
    groups = ["stress only", "brine", "melt"]
    fig = Figure(size = (1350, 560), fontsize = 14)
    Label(fig[0, 1:3], "Angle between turbulent stress and shear in the LES (0° = down-gradient, as every column scheme assumes)";
          fontsize = 16, font = :bold)
    for (p, g) in enumerate(groups)
        ax = Axis(fig[1, p]; title = g, xlabel = "θ (degrees, > 0: stress to the left of the shear)",
                  ylabel = p == 1 ? "depth / mixed-layer depth" : "", yreversed = true)
        vlines!(ax, [0]; color = (:black, 0.4))
        hlines!(ax, [1]; color = (:black, 0.3), linestyle = :dot)
        members = [r for r in results if r.flux == g]
        colors = Makie.wong_colors()
        for (i, r) in enumerate(members)
            ζ = [a.ζ for a in r.angle_profile if a.n ≥ 5]
            med = [a.median for a in r.angle_profile if a.n ≥ 5]
            isempty(ζ) && continue
            lines!(ax, med, ζ; color = colors[mod1(i, 7)], linewidth = 2.5,
                   linestyle = r.f == 0 ? :dash : :solid, label = replace(r.name, "arctic_ice_drag" => "ice drag"))
        end
        isempty(members) || axislegend(ax; position = :rb, labelsize = 10, framevisible = false)
        xlims!(ax, -90, 90)
    end
    save(joinpath(folder, "stress_shear_angle.png"), fig; px_per_unit = 2)

    # Figure 2: h f / u★ against inertial periods (rotating cases), and h against h_PRT
    fig2 = Figure(size = (1250, 520), fontsize = 14)
    ax1 = Axis(fig2[1, 1]; xlabel = "inertial periods, f t / 2π", ylabel = "h |f| / u★",
               title = "Rotating cases: neutral Ekman scaling predicts 0.25–0.4")
    band!(ax1, [0, 10], [0.25, 0.25], [0.4, 0.4]; color = (:grey, 0.2))
    ax2 = Axis(fig2[1, 2]; xlabel = "1.7 u★ / √(N f)  (m)", ylabel = "LES h (m)",
               title = "Pollard–Rhines–Thompson (second-half medians)")
    colors = Makie.wong_colors()
    markers = Dict("stress only" => :circle, "brine" => :utriangle, "melt" => :dtriangle)
    rotating = [r for r in results if r.f != 0]
    for (i, r) in enumerate(rotating)
        c = colors[mod1(i, 7)]
        lines!(ax1, abs(r.f) .* r.t ./ 2π, r.h_t .* abs(r.f) ./ max.(r.u★_t, 1e-4); color = c, linewidth = 2,
               label = replace(r.name, "arctic_ice_drag" => "ice drag"))
        scatter!(ax2, [r.prt], [r.h]; color = c, marker = markers[r.flux], markersize = 14)
    end
    upper = maximum(vcat([r.prt for r in rotating], [r.h for r in rotating], [10.0]); init = 10.0) * 1.1
    lines!(ax2, [0, upper], [0, upper]; color = (:black, 0.4), linestyle = :dash)
    isempty(rotating) || axislegend(ax1; position = :rt, labelsize = 9, framevisible = false)
    ylims!(ax1, 0, 1.2)
    save(joinpath(folder, "depth_scaling.png"), fig2; px_per_unit = 2)
    return results
end

rotation_results = run_rotation_diagnostics(output_folder())
@info "Wrote $(joinpath(output_folder(), "rotation_diagnostics.md")), stress_shear_angle.png and depth_scaling.png"
