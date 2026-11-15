# Experimental fixes for the bottom-cell pressure gradient on z-star grids, as run-time method overrides.
# Nothing in the installed Oceananigans is changed: `include` this file after `using Oceananigans` and the new
# methods exist for this Julia session only. Each fix is off until its switch is set, so one session can run every
# case; with both switches off the model is exactly the stock one.
#
# Written against Oceananigans fjp/net-inflow-target (53efe85, v0.117.2 + NetInflow):
#   src/ImmersedBoundaries/mutable_immersed_grid.jl             z-star chain rule  ∂xᶠᶜᶜ, ∂yᶜᶠᶜ (grid::AMG, ϕ)
#   src/Models/NonhydrostaticModels/update_hydrostatic_pressure.jl   pHY′ integration
#
# What the stock chain-rule term sees in a bottom cell k = kᵇ (∂z p′ averaged over faces kᵇ and kᵇ+1 of both columns):
#   face kᵇ+1: δz p′ / (partial spacing) = ℑz b × (full / partial spacing), up to ×1.7 for φ = 0.2;
#   face kᵇ:   0. On the immersed grid δz is masked at faces next to an immersed cell (conditional_differences.jl),
#              so the bottom face contributes nothing instead of b. (Found 15 Nov with probe_fixes.jl: F2 alone
#              changes p′ in immersed cells by up to 624 m²/s² but the force not at all.)
#
# F1 — chain-rule ∂z p′ taken on the underlying grid: full spacing and no masking. Face kᵇ+1 then gives ℑz b exactly.
#   Face kᵇ gives δz p′ between the bottom cell and the immersed cell below it, i.e. ½(b_kᵇ + b_immersed): fresh
#   0 °C water unless F2 is also on. So F1 alone replaces 0 by ½(b + b₀) at the bottom face.
#   Partial- and whole-cell grids (whole cells have the bottom-face masking too).
#
# F2 — immersed cells take the bottom wet cell's buoyancy in the p′ integration. Wet cells' p′ is unchanged; only p′
#   in immersed cells changes, which matters only through F1's unmasked bottom face. F1+F2: ∂z p′ = b at both faces.
#
# Inside a kernel every Field arrives as an offset array (`adapt`, also on the CPU; on 0.117 a
# StaticStridesOffsetArray), so F1 is defined for those. In the hydrostatic model the only ∂xᶠᶜᶜ/∂yᶜᶠᶜ of an array is
# the pressure gradient (other callers are closures and the implicit free surface, not used in the ridge test).
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

# Kernel arguments: on 0.117 a Field arrives in a kernel as a StaticStridesOffsetArray (plain OffsetArray on others).
const OffsetArray3 = Union{Oceananigans.Grids.OffsetArrays.OffsetArray{<:Any, 3},
                           Oceananigans.Utils.StaticStridesOffsetArray{<:Any, 3}}

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

const BottomMutableGrid = Union{PCBMutableGrid, GFBMutableGrid}

@inline function Operators.∂xᶠᶜᶜ(i, j, k, grid::BottomMutableGrid, ϕ::OffsetArray3)
    ∂x_at_r = δxᶠᶜᶜ(i, j, k, grid, ϕ) * Δx⁻¹ᶠᶜᶜ(i, j, k, grid)
    if FIX_F1[]
        Threads.atomic_add!(F1_CALLS, 1)
        ∂z_ϕ = ℑxzᶠᵃᶜ(i, j, k, grid.underlying_grid, ∂zᶜᶜᶠ, ϕ)
    else
        ∂z_ϕ = ℑxzᶠᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ϕ)
    end
    return ∂x_at_r - ∂x_zᶠᶜᶜ(i, j, k, grid) * ∂z_ϕ
end

@inline function Operators.∂yᶜᶠᶜ(i, j, k, grid::BottomMutableGrid, ϕ::OffsetArray3)
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
