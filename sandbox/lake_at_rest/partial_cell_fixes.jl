# Experimental fixes for the bottom-cell pressure gradient on z-star grids, as run-time method overrides.
# Nothing in the installed Oceananigans is changed: `include` this file after `using Oceananigans` and the new
# methods exist for this Julia session only. Each fix is off until its switch is set, so one session can run every
# case; with both switches off the model is exactly the stock one.
#
# Written against Oceananigans fjp/net-inflow-target (53efe85, v0.117.2 + NetInflow):
#   src/ImmersedBoundaries/mutable_immersed_grid.jl             z-star chain rule  ∂xᶠᶜᶜ, ∂yᶜᶠᶜ (grid::AMG, ϕ)
#   src/Models/NonhydrostaticModels/update_hydrostatic_pressure.jl   pHY′ integration
#
# F1 — chain-rule ∂z p′ with the full-cell face spacing.
#   The force is −∂x p′|_z = −δx p′/Δx + (∂x z)(∂z p′), with ∂z p′ = ℑxz(∂zᶜᶜᶠ p′). On a partial-cell grid ∂zᶜᶜᶠ divides
#   by the partial face spacing (wet midpoint, up to ×1.7 at face kᵇ+1), while p′ was integrated with the full
#   spacing. F1 divides by the full spacing (the underlying grid's), so that ∂z p′ = ℑz b at every wet face, as the
#   integration says. Partial-cell grids only (on whole-cell grids the two spacings are the same).
#
# F2 — immersed cells take the bottom wet cell's buoyancy in the p′ integration.
#   Immersed cells hold T = S = 0, and the integration runs through them, so p′ in the first immersed cell (used by
#   ∂z p′ at face kᵇ in the chain rule) is built from fresh 0 °C water. F2 extrapolates the last wet cell's buoyancy
#   downward instead. Wet cells are unchanged. Both partial- and whole-cell grids (the integral passes through
#   immersed cells in both).
#
# Inside a kernel every Field arrives as its OffsetArray (`adapt`, also on the CPU), so F1 is defined for 3-D
# OffsetArrays. In the hydrostatic model the only ∂xᶠᶜᶜ/∂yᶜᶠᶜ of an array is the pressure gradient (the other callers
# are closures and the implicit free surface, not used in the ridge test), so F1 touches nothing else there.
# `F1_CALLS[]` counts the overridden calls: it must be > 0 after a step with F1 on.

using Oceananigans
using Oceananigans.Grids: Center, Bounded, inactive_node
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, PartialCellVerticalBoundaries,
                                       GridFittedVerticalBoundaries, UnderlyingMutableGrid
using Oceananigans.Operators: Operators, δxᶠᶜᶜ, δyᶜᶠᶜ, Δx⁻¹ᶠᶜᶜ, Δy⁻¹ᶜᶠᶜ, ℑxzᶠᵃᶜ, ℑyzᵃᶠᶜ, ∂zᶜᶜᶠ,
                               ∂x_zᶠᶜᶜ, ∂y_zᶜᶠᶜ, Δzᶜᶜᶜ, Δzᶜᶜᶠ
using Oceananigans.BuoyancyFormulations: z_dot_g_bᶜᶜᶜ
using Oceananigans.Utils: launch!
using Oceananigans.Grids: surface_kernel_parameters
using Oceananigans.Models.NonhydrostaticModels: @kernel, @index   # KernelAbstractions, via Oceananigans (not a direct dependency of the nested env)

import Oceananigans.Models.NonhydrostaticModels: update_hydrostatic_pressure!

const OffsetArray3 = Oceananigans.Grids.OffsetArrays.OffsetArray{<:Any, 3}

const FIX_F1 = Ref(false)
const FIX_F2 = Ref(false)
const F1_CALLS = Threads.Atomic{Int}(0)

function set_fixes!(; F1 = false, F2 = false)
    FIX_F1[] = F1
    FIX_F2[] = F2
    F1_CALLS[] = 0
    return nothing
end

fix_label() = FIX_F1[] & FIX_F2[] ? "F1+F2" : FIX_F1[] ? "F1" : FIX_F2[] ? "F2" : "none"

# Mutable (z-star) immersed grids with a partial-cell or whole-cell bottom. More specific than every stock method
# they shadow, so no stock method is overwritten.
const PCBMutableGrid = ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Bounded, <:UnderlyingMutableGrid, <:PartialCellVerticalBoundaries}
const GFBMutableGrid = ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Bounded, <:UnderlyingMutableGrid, <:GridFittedVerticalBoundaries}

#####
##### F1
#####

@inline function Operators.∂xᶠᶜᶜ(i, j, k, grid::PCBMutableGrid, ϕ::OffsetArray3)
    ∂x_at_r = δxᶠᶜᶜ(i, j, k, grid, ϕ) * Δx⁻¹ᶠᶜᶜ(i, j, k, grid)
    if FIX_F1[]
        Threads.atomic_add!(F1_CALLS, 1)
        ∂z_ϕ = ℑxzᶠᵃᶜ(i, j, k, grid.underlying_grid, ∂zᶜᶜᶠ, ϕ)
    else
        ∂z_ϕ = ℑxzᶠᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ϕ)
    end
    return ∂x_at_r - ∂x_zᶠᶜᶜ(i, j, k, grid) * ∂z_ϕ
end

@inline function Operators.∂yᶜᶠᶜ(i, j, k, grid::PCBMutableGrid, ϕ::OffsetArray3)
    ∂y_at_r = δyᶜᶠᶜ(i, j, k, grid, ϕ) * Δy⁻¹ᶜᶠᶜ(i, j, k, grid)
    if FIX_F1[]
        Threads.atomic_add!(F1_CALLS, 1)
        ∂z_ϕ = ℑyzᵃᶠᶜ(i, j, k, grid.underlying_grid, ∂zᶜᶜᶠ, ϕ)
    else
        ∂z_ϕ = ℑyzᵃᶠᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ϕ)
    end
    return ∂y_at_r - ∂y_zᶜᶠᶜ(i, j, k, grid) * ∂z_ϕ
end

#####
##### F2
#####

# Same integration as the stock kernel (half a top cell, then the trapezoidal rule with full spacings), except that an
# immersed cell takes the buoyancy of the cell above it, so the last wet cell's buoyancy continues downward.
@kernel function _hydrostatic_pressure_filled!(pHY′, ibg, buoyancy, C)
    i, j = @index(Global, NTuple)
    grid = ibg.underlying_grid
    Nz = grid.Nz

    b_above = z_dot_g_bᶜᶜᶜ(i, j, Nz, grid, buoyancy, C)
    pᵏ = - b_above * Δzᶜᶜᶜ(i, j, Nz, grid) / 2
    @inbounds pHY′[i, j, Nz] = pᵏ

    for k in Nz-1 : -1 : 1
        wet = !inactive_node(i, j, k, ibg, Center(), Center(), Center())
        b_here = ifelse(wet, z_dot_g_bᶜᶜᶜ(i, j, k, grid, buoyancy, C), b_above)
        pᵏ -= (b_here + b_above) / 2 * Δzᶜᶜᶠ(i, j, k+1, grid)
        @inbounds pHY′[i, j, k] = pᵏ
        b_above = b_here
    end
end

# Stock behaviour when the switch is off: partial-cell grids integrate on the underlying grid, whole-cell grids on
# the immersed grid itself (see update_hydrostatic_pressure.jl).
function update_hydrostatic_pressure!(pHY′, arch, ibg::PCBMutableGrid, buoyancy, tracers;
                                      parameters = surface_kernel_parameters(ibg.underlying_grid))
    if FIX_F2[]
        launch!(arch, ibg.underlying_grid, parameters, _hydrostatic_pressure_filled!, pHY′, ibg, buoyancy, tracers)
    else
        update_hydrostatic_pressure!(pHY′, arch, ibg.underlying_grid, buoyancy, tracers; parameters)
    end
    return nothing
end

function update_hydrostatic_pressure!(pHY′, arch, ibg::GFBMutableGrid, buoyancy, tracers;
                                      parameters = surface_kernel_parameters(ibg))
    if FIX_F2[]
        launch!(arch, ibg.underlying_grid, parameters, _hydrostatic_pressure_filled!, pHY′, ibg, buoyancy, tracers)
    else
        update_hydrostatic_pressure!(pHY′, arch, ibg.underlying_grid, buoyancy, tracers; parameters)
    end
    return nothing
end
