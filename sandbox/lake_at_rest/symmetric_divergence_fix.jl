# Candidate fix for WENOVectorInvariant at partial cells, as a run-time override (Oceananigans itself is unchanged).
#
# WENOVectorInvariant writes the vertical advection as
#     w ∂z u = V⁻¹ [ δz(W u) + u (horizontal divergence at the u-face) ]
# and takes that divergence by WENO-biased (upwind) interpolation of the neighbouring T-cells' divergences
# (upwinded_divergence_flux_Uᶠᶜᶜ in Advection/vector_invariant_self_upwinding.jl). At a thin partial cell next to a
# thick one, the upwind cell can be the thin one with almost no divergence while the thick one carries most of the
# u-cell's mass budget, so the two parts no longer cancel and a bottom jet forms (mwe_partial_cells.jl).
#
# The fix: interpolate the divergence symmetrically (the mean of the two neighbouring T-cells, which is the u-cell's
# own horizontal divergence δx ℑx(Ax u)). Gated by SYMMETRIC_DIVERGENCE[] so that one session can run both.
# Applied to every immersed grid (whole and partial cells) so the comparison stays fair.

using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using Oceananigans.Advection: Advection, VectorInvariantSelfVerticalUpwinding, upwinded_divergence_flux_Uᶠᶜᶜ,
                              upwinded_divergence_flux_Vᶜᶠᶜ, vector_invariant_flux_Wu, vector_invariant_flux_Wv,
                              δx_U_plus_∂t_σ, δy_V
using Oceananigans.Operators: V⁻¹ᶠᶜᶜ, V⁻¹ᶜᶠᶜ, δzᵃᵃᶜ, ℑxᶠᵃᵃ, ℑyᵃᶠᵃ

const SYMMETRIC_DIVERGENCE = Ref(false)

# Horizontal divergence of a T-cell (plus the z-star grid-motion term, as in the stock self-upwinding)
@inline cell_divergence(i, j, k, grid, u, v) = δx_U_plus_∂t_σ(i, j, k, grid, u, v) + δy_V(i, j, k, grid, u, v)

@inline function Advection.vertical_advection_U(i, j, k, grid::ImmersedBoundaryGrid,
                                                scheme::VectorInvariantSelfVerticalUpwinding, U)
    𝒜ᶻ = δzᵃᵃᶜ(i, j, k, grid, vector_invariant_flux_Wu, scheme, U.w, U.u)
    Φᵟ = SYMMETRIC_DIVERGENCE[] ?
         (@inbounds U.u[i, j, k]) * ℑxᶠᵃᵃ(i, j, k, grid, cell_divergence, U.u, U.v) :
         upwinded_divergence_flux_Uᶠᶜᶜ(i, j, k, grid, scheme, U.u, U.v)
    return V⁻¹ᶠᶜᶜ(i, j, k, grid) * (Φᵟ + 𝒜ᶻ)
end

@inline function Advection.vertical_advection_V(i, j, k, grid::ImmersedBoundaryGrid,
                                                scheme::VectorInvariantSelfVerticalUpwinding, U)
    𝒜ᶻ = δzᵃᵃᶜ(i, j, k, grid, vector_invariant_flux_Wv, scheme, U.w, U.v)
    Φᵟ = SYMMETRIC_DIVERGENCE[] ?
         (@inbounds U.v[i, j, k]) * ℑyᵃᶠᵃ(i, j, k, grid, cell_divergence, U.u, U.v) :
         upwinded_divergence_flux_Vᶜᶠᶜ(i, j, k, grid, scheme, U.u, U.v)
    return V⁻¹ᶜᶠᶜ(i, j, k, grid) * (Φᵟ + 𝒜ᶻ)
end
