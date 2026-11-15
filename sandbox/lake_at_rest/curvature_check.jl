# Is the hybrid fix's missing curvature negligible on the pan-Arctic grid?
#
# At faces next to a bottom-thickness change the hybrid fix uses the flux form of momentum advection, which in
# Oceananigans has no curvature (metric) terms on curvilinear grids; the vector-invariant form includes them
# implicitly. This script estimates, at every face where the fix acts, the size of the missing terms from the model's
# current state and compares it with the Coriolis term there. For an orthogonal grid with spacings Δx, Δy (per grid
# step), the curvature terms of the momentum equations are
#
#     u:  v̄ (v̄ δx Δy − u δy Δx) / (Δx Δy)        v:  ū (ū δy Δx − v δx Δy) / (Δx Δy)
#
# (zero on a rectilinear grid). The script also gives the difference between the hybrid's and the stock scheme's
# advection at the same faces, i.e. the change the fix makes on purpose.
#
# Use in the REPL of a finished hybrid run (the session that ran set_momentum_fix!("hybrid")), with the ocean's
# model, e.g.
#
#     include("<panarctic-tools>/sandbox/lake_at_rest/curvature_check.jl")
#     curvature_check(ocean.model)      # whatever panarctic_nested.jl calls the ocean model
#
# CPU loop over a copy of the state; a minute or two at 15 km.

isdefined(Main, :thickness_changes_x) || include(joinpath(@__DIR__, "thickness_jump_fixes.jl"))

using Printf
using Oceananigans.Architectures: on_architecture, CPU
using Oceananigans.Grids: peripheral_node, Center, Face
using Oceananigans.Operators: Δxᶠᶜᶜ, Δyᶠᶜᶜ, Δxᶜᶠᶜ, Δyᶜᶠᶜ, Δxᶜᶜᶜ, Δyᶜᶜᶜ, Δxᶠᶠᶜ, Δyᶠᶠᶜ, ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using Oceananigans.Coriolis: x_f_cross_U, y_f_cross_U
using Oceananigans.Advection: Advection

@inline function curvature_u(i, j, k, grid, u, v)
    v̄ = ℑxyᶠᶜᵃ(i, j, k, grid, v)
    δxΔy = Δyᶜᶜᶜ(i, j, k, grid) - Δyᶜᶜᶜ(i - 1, j, k, grid)
    δyΔx = Δxᶠᶠᶜ(i, j + 1, k, grid) - Δxᶠᶠᶜ(i, j, k, grid)
    return v̄ * (v̄ * δxΔy - (@inbounds u[i, j, k]) * δyΔx) / (Δxᶠᶜᶜ(i, j, k, grid) * Δyᶠᶜᶜ(i, j, k, grid))
end

@inline function curvature_v(i, j, k, grid, u, v)
    ū = ℑxyᶜᶠᵃ(i, j, k, grid, u)
    δyΔx = Δxᶜᶜᶜ(i, j, k, grid) - Δxᶜᶜᶜ(i, j - 1, k, grid)
    δxΔy = Δyᶠᶠᶜ(i + 1, j, k, grid) - Δyᶠᶠᶜ(i, j, k, grid)
    return ū * (ū * δyΔx - (@inbounds v[i, j, k]) * δxΔy) / (Δxᶜᶠᶜ(i, j, k, grid) * Δyᶜᶠᶜ(i, j, k, grid))
end

# Stock vector-invariant advection (through the override code; with the "hybrid" switches these are stock)
stock_u(i, j, k, grid, s, U) = Advection.horizontal_advection_U(i, j, k, grid, s, U.u, U.v) +
                               Advection.vertical_advection_U(i, j, k, grid, s, U) +
                               Advection.bernoulli_head_U(i, j, k, grid, s, U.u, U.v)
stock_v(i, j, k, grid, s, U) = Advection.horizontal_advection_V(i, j, k, grid, s, U.u, U.v) +
                               Advection.vertical_advection_V(i, j, k, grid, s, U) +
                               Advection.bernoulli_head_V(i, j, k, grid, s, U.u, U.v)

quantile_of(x, p) = isempty(x) ? NaN : sort(x)[clamp(ceil(Int, p * length(x)), 1, length(x))]

function curvature_check(model; top = 10)
    (ke_fix_on() || symmetric_divergence_on()) &&
        @warn "ke/symdiv switches are on: the 'stock' column below is not stock. Run in the hybrid session."
    grid = on_architecture(CPU(), model.grid)
    U = map(f -> on_architecture(CPU(), f), model.velocities)
    coriolis = on_architecture(CPU(), model.coriolis)
    scheme = model.advection.momentum
    flux = scheme.kinetic_energy_gradient_scheme
    Nx, Ny, Nz = size(grid)
    C, F = Center(), Face()
    rows = NamedTuple[]
    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        if !peripheral_node(i, j, k, grid, F, C, C) && thickness_changes_x(i, j, k, grid)
            M = curvature_u(i, j, k, grid, U.u, U.v)
            f = x_f_cross_U(i, j, k, grid, coriolis, U)
            D = Advection.div_𝐯u(i, j, k, grid, flux, U, U.u) - stock_u(i, j, k, grid, scheme, U)
            push!(rows, (; face = :u, i, j, k, M = abs(M), C = abs(f), D = abs(D), speed = abs(U.u[i, j, k])))
        end
        if !peripheral_node(i, j, k, grid, C, F, C) && thickness_changes_y(i, j, k, grid)
            M = curvature_v(i, j, k, grid, U.u, U.v)
            f = y_f_cross_U(i, j, k, grid, coriolis, U)
            D = Advection.div_𝐯v(i, j, k, grid, flux, U, U.v) - stock_v(i, j, k, grid, scheme, U)
            push!(rows, (; face = :v, i, j, k, M = abs(M), C = abs(f), D = abs(D), speed = abs(U.v[i, j, k])))
        end
    end
    n = length(rows)
    ratio = [r.M / max(r.C, 1e-12) for r in rows]
    Dratio = [r.M / max(r.D, 1e-15) for r in rows]
    perday = 86400 * 100                                    # m/s² → cm/s per day
    @printf("Faces where the hybrid fix acts: %d (walls %s)\n\n", n, jump_counts_walls() ? "counted" : "not counted")
    println("Missing curvature term M, against the Coriolis term C and the hybrid − stock difference D, at those faces:")
    @printf("  M:     median %.2e, 99%% %.2e, max %.2e m/s²  (max = %.2f cm/s per day)\n",
            quantile_of([r.M for r in rows], 0.5), quantile_of([r.M for r in rows], 0.99), maximum(r.M for r in rows),
            perday * maximum(r.M for r in rows))
    @printf("  M / C: median %.2e, 99%% %.2e, max %.2e;  faces with M > 1 %% of C: %.2f %%\n",
            quantile_of(ratio, 0.5), quantile_of(ratio, 0.99), maximum(ratio), 100count(>(0.01), ratio) / n)
    @printf("  M / D: median %.2e, 99%% %.2e  (D median %.2e m/s²: the change the fix makes on purpose)\n",
            quantile_of(Dratio, 0.5), quantile_of(Dratio, 0.99), quantile_of([r.D for r in rows], 0.5))
    println("\nLargest M / C (where it matters most):")
    @printf("  %-4s %5s %5s %3s  %-10s %-10s %-10s %-9s  %s\n", "face", "i", "j", "k", "M", "C", "D", "speed", "lon, lat")
    for r in sort(rows, by = r -> -r.M / max(r.C, 1e-12))[1:min(top, n)]
        ℓ = r.face == :u ? (F, C, C) : (C, F, C)
        place = try
            @sprintf("%.1f°E %.1f°N", Oceananigans.Grids.λnode(r.i, r.j, r.k, grid, ℓ...),
                     Oceananigans.Grids.φnode(r.i, r.j, r.k, grid, ℓ...))
        catch
            ""
        end
        @printf("  %-4s %5d %5d %3d  %-10.2e %-10.2e %-10.2e %-9.3f  %s\n", r.face, r.i, r.j, r.k, r.M, r.C, r.D, r.speed, place)
    end
    println("""

Reading it: if M is below ~1 % of C at nearly all faces and well below D, the missing curvature is negligible next to
both the rotation and the change the fix makes, and the hybrid is fine on this grid. If the largest M / C cluster
in one place (near the projection's pole, or where the grid spacing changes fastest), that is where a vector-invariant
fix ("ke") would be safer.""")
    return rows
end
