# Sandbox: M5b ("balanced Flather") for the pan-Arctic model, written against panarctic.jl as of 3 Nov 2026
# (Oceananigans 0.113.6/0.116.1). The model has since moved its setup into setup.jl; check the hook points
# (atlantic_V_obc, pacific_V_obc, V_obcs) before using it. See sandbox/obc-vortex/README.md for the tests.

# Balanced Flather (M5b) for panarctic.jl: Flather about ηᵉˣᵗ + C on each open edge, with C chosen every barotropic
# substep so that the net transport through that edge equals the external (GLORYS) net transport:
#
#     V = Vᵉˣᵗ ∓ c (ηᵇ − ηᵉˣᵗ − C),    C = Σ c (ηᵇ − ηᵉˣᵗ) Δx / Σ c Δx,    c = √(g H)
#
# so a uniform offset between the model's and GLORYS's free surface along an edge drives no transport, and the
# correction is spread as √(gH) (gentle on shelves), unlike target_transport's uniform shift per unit length.
# Tested in obc-vortex (claude/obc-vortex-results.md).
#
# Used by panarctic_setup_bf.jl, which includes this file just before V_obcs is built and swaps atlantic_V_obc /
# pacific_V_obc for the versions below. Needs atlantic_V_obc, pacific_V_obc and ηᵉˣᵗ from the script.

using Oceananigans.Grids: static_column_depthᶜᶠᵃ
using Oceananigans.Operators: Δxᶜᶠᶜ

const g_balanced_flather = Oceananigans.defaults.gravitational_acceleration

# C for one edge. jᶠ: the boundary face row; jᵢ, jₕ: the interior and halo rows of η either side of it;
# jₑ: the row of ηᵉˣᵗ the script compares with (1 at the south, Ny at the north).
@inline function edge_offset(jᶠ, jᵢ, jₕ, jₑ, grid, clock, fields, ηᵉˣᵗ)
    t  = Oceananigans.Units.Time(clock.time)
    η  = fields.η
    kη = grid.Nz + 1
    num = zero(grid)
    den = zero(grid)
    for i′ in 1:grid.Nx
        H  = static_column_depthᶜᶠᵃ(i′, jᶠ, grid)
        c  = sqrt(g_balanced_flather * max(H, zero(H)))
        Δx = Δxᶜᶠᶜ(i′, jᶠ, 1, grid)
        ηᵇ = (η[i′, jᵢ, kη] + η[i′, jₕ, kη]) / 2
        num += ifelse(H > 0, c * (ηᵇ - ηᵉˣᵗ[i′, jₑ, 1, t]) * Δx, zero(grid))
        den += c * Δx
    end
    return ifelse(den > 0, num / den, zero(grid))
end

@inline function atlantic_V_obc_balanced(i, k, grid, clock, fields, p)
    V, ηᵉ = atlantic_V_obc(i, k, grid, clock, fields, p)
    C = edge_offset(1, 1, 0, 1, grid, clock, fields, p.η)
    return (V, ηᵉ + C)
end

@inline function pacific_V_obc_balanced(i, k, grid, clock, fields, p)
    V, ηᵉ = pacific_V_obc(i, k, grid, clock, fields, p)
    C = edge_offset(grid.Ny + 1, grid.Ny, grid.Ny + 1, grid.Ny, grid, clock, fields, p.η)
    return (V, ηᵉ + C)
end

@info "Balanced Flather (M5b) on the Atlantic and Bering edges"
