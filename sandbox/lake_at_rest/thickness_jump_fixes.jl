# Two candidate fixes for WENOVectorInvariant at bottom cells whose thickness changes from one column to the next
# (thin partial cells next to thick ones, and cliff edges). Run-time overrides: Oceananigans itself is unchanged.
#
# Background (handoff-2026-11-15b): the flux form of momentum advection brakes the jet that forms in a thin bottom
# cell at a cliff, because each u-cell's momentum budget uses the same transports as its mass budget. The
# vector-invariant form splits advection into a kinetic-energy gradient, a vertical term and a vorticity term, and
# two of those lose the u-cell's mass budget where the cell thickness jumps:
#   (a) the vertical term uses an upwinded divergence (fixed by symmetric_divergence_fix.jl, "symdiv");
#   (b) the kinetic-energy gradient ∂x(u²/2) treats a thin and a thick face as equals.
#
# "ke" fix (on top of symdiv): at faces where the thickness changes, replace ∂x(u²/2) by the advective term the flux
# form implies,
#     u ∂x u  ≈  V⁻¹ [ δx(Ū ũ) − u δx Ū ],        Ū = ℑx(Ax u) (the transport through each T-cell centre),
# where Ū ũ is the flux form's momentum flux (WENO, the scheme's own kinetic-energy-gradient scheme). With symdiv's
# divergence (the u-cell's own, ℑx of the two T-cell divergences) the two parts add up to the flux form at those faces.
#
# "hybrid" fix: at faces where the thickness changes, use the flux form for the whole advection term; elsewhere the
# stock vector-invariant form. Nothing else changes. (At these faces the flux form uses no curvature terms: fine on
# the rectilinear test grid; to be checked on the pan-Arctic grid.)
#
# "Where the thickness changes": the u-face's stencil reaches the four T-cells i−2 … i+1 (the transports at the two
# centres i−1 and i use faces i−1, i, i+1, and each face's area uses its two cells). The fix acts if their wet
# thicknesses differ by more than 1 % (an immersed cell counts as 0, so cliff edges are included). In the interior of
# a flat-bottomed or smoothly varying region nothing changes.
#
# Switches are functions, not Refs, so the same file works in GPU kernels. Use set_momentum_fix!(name) and run the
# model through Base.invokelatest (the switch redefines a method).

include(joinpath(@__DIR__, "symmetric_divergence_fix.jl"))   # symdiv, used by "ke"

using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using Oceananigans.Grids: inactive_node, topology, Flat, Center, Face
using Oceananigans.Advection: Advection, VectorInvariantKineticEnergyUpwinding,
                              _advective_momentum_flux_Uu, _advective_momentum_flux_Vv,
                              δx_u², δx_v², δy_u², δy_v², bias,
                              _biased_interpolate_xᶠᵃᵃ, _biased_interpolate_yᵃᶠᵃ,
                              _symmetric_interpolate_yᵃᶜᵃ, _symmetric_interpolate_xᶜᵃᵃ
using Oceananigans.Operators: V⁻¹ᶠᶜᶜ, V⁻¹ᶜᶠᶜ, Δx⁻¹ᶠᶜᶜ, Δy⁻¹ᶜᶠᶜ, δxᶠᵃᵃ, δyᵃᶠᵃ, ℑxᶜᵃᵃ, ℑyᵃᶜᵃ,
                              Ax_qᶠᶜᶜ, Ay_qᶜᶠᶜ, Δrᶜᶜᶜ

# Switches (redefined by set_momentum_fix!)
ke_fix_on() = false
hybrid_on() = false

#####
##### Where the thickness changes
#####

@inline wet_thickness(i, j, k, grid) =
    ifelse(inactive_node(i, j, k, grid, Center(), Center(), Center()), zero(grid), Δrᶜᶜᶜ(i, j, k, grid))

@inline function thickness_changes(a, b, c, d)
    hi = max(a, b, c, d)
    lo = min(a, b, c, d)
    return hi - lo > hi / 100
end

@inline thickness_changes_x(i, j, k, grid) =
    thickness_changes(wet_thickness(i - 2, j, k, grid), wet_thickness(i - 1, j, k, grid),
                      wet_thickness(i,     j, k, grid), wet_thickness(i + 1, j, k, grid))

@inline thickness_changes_y(i, j, k, grid) =
    thickness_changes(wet_thickness(i, j - 2, k, grid), wet_thickness(i, j - 1, k, grid),
                      wet_thickness(i, j,     k, grid), wet_thickness(i, j + 1, k, grid))

# No y-direction fix on grids that are Flat in y (the 2-D test): there are no neighbours in y.
@inline acts_in_y(grid) = !(topology(grid, 2) === Flat)

#####
##### "ke": the flux form's u ∂x u (and v ∂y v) at faces where the thickness changes
#####

@inline Ūᶜᶜᶜ(i, j, k, grid, u) = ℑxᶜᵃᵃ(i, j, k, grid, Ax_qᶠᶜᶜ, u)
@inline V̄ᶜᶜᶜ(i, j, k, grid, v) = ℑyᵃᶜᵃ(i, j, k, grid, Ay_qᶜᶠᶜ, v)

@inline function flux_form_u∂xu(i, j, k, grid, weno, u)
    F = δxᶠᵃᵃ(i, j, k, grid, _advective_momentum_flux_Uu, weno, u, u)   # δx(Ū ũ)
    D = δxᶠᵃᵃ(i, j, k, grid, Ūᶜᶜᶜ, u)                                  # δx Ū
    return V⁻¹ᶠᶜᶜ(i, j, k, grid) * (F - (@inbounds u[i, j, k]) * D)
end

@inline function flux_form_v∂yv(i, j, k, grid, weno, v)
    F = δyᵃᶠᵃ(i, j, k, grid, _advective_momentum_flux_Vv, weno, v, v)
    D = δyᵃᶠᵃ(i, j, k, grid, V̄ᶜᶜᶜ, v)
    return V⁻¹ᶜᶠᶜ(i, j, k, grid) * (F - (@inbounds v[i, j, k]) * D)
end

# Same as the stock bernoulli_head_U for WENOVectorInvariant, except at thickness changes when "ke" is on
@inline function Advection.bernoulli_head_U(i, j, k, grid::ImmersedBoundaryGrid,
                                            scheme::VectorInvariantKineticEnergyUpwinding, u, v)
    @inbounds û = u[i, j, k]
    cross_scheme = scheme.upwinding.cross_scheme
    δKvˢ = _symmetric_interpolate_yᵃᶜᵃ(i, j, k, grid, scheme, cross_scheme, δx_v², u, v)   # v ∂x v part: unchanged
    if ke_fix_on() && thickness_changes_x(i, j, k, grid)
        return δKvˢ * Δx⁻¹ᶠᶜᶜ(i, j, k, grid) + flux_form_u∂xu(i, j, k, grid, scheme.kinetic_energy_gradient_scheme, u)
    else
        δKuᴿ = _biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, scheme.kinetic_energy_gradient_scheme, bias(û),
                                        δx_u², scheme.upwinding.δu²_stencil, u, v)
        return (δKuᴿ + δKvˢ) * Δx⁻¹ᶠᶜᶜ(i, j, k, grid)
    end
end

@inline function Advection.bernoulli_head_V(i, j, k, grid::ImmersedBoundaryGrid,
                                            scheme::VectorInvariantKineticEnergyUpwinding, u, v)
    @inbounds v̂ = v[i, j, k]
    cross_scheme = scheme.upwinding.cross_scheme
    δKuˢ = _symmetric_interpolate_xᶜᵃᵃ(i, j, k, grid, scheme, cross_scheme, δy_u², u, v)
    if ke_fix_on() && acts_in_y(grid) && thickness_changes_y(i, j, k, grid)
        return δKuˢ * Δy⁻¹ᶜᶠᶜ(i, j, k, grid) + flux_form_v∂yv(i, j, k, grid, scheme.kinetic_energy_gradient_scheme, v)
    else
        δKvᴿ = _biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, scheme.kinetic_energy_gradient_scheme, bias(v̂),
                                        δy_v², scheme.upwinding.δv²_stencil, u, v)
        return (δKvᴿ + δKuˢ) * Δy⁻¹ᶜᶠᶜ(i, j, k, grid)
    end
end

#####
##### "hybrid": the flux form for the whole advection term at faces where the thickness changes
#####

@inline function Advection.U_dot_∇u(i, j, k, grid::ImmersedBoundaryGrid,
                                    scheme::VectorInvariantKineticEnergyUpwinding, U)
    if hybrid_on() && thickness_changes_x(i, j, k, grid)
        return Advection.div_𝐯u(i, j, k, grid, scheme.kinetic_energy_gradient_scheme, U, U.u)
    else
        return Advection.horizontal_advection_U(i, j, k, grid, scheme, U.u, U.v) +
               Advection.vertical_advection_U(i, j, k, grid, scheme, U) +
               Advection.bernoulli_head_U(i, j, k, grid, scheme, U.u, U.v)
    end
end

@inline function Advection.U_dot_∇v(i, j, k, grid::ImmersedBoundaryGrid,
                                    scheme::VectorInvariantKineticEnergyUpwinding, U)
    if hybrid_on() && acts_in_y(grid) && thickness_changes_y(i, j, k, grid)
        return Advection.div_𝐯v(i, j, k, grid, scheme.kinetic_energy_gradient_scheme, U, U.v)
    else
        return Advection.horizontal_advection_V(i, j, k, grid, scheme, U.u, U.v) +
               Advection.vertical_advection_V(i, j, k, grid, scheme, U) +
               Advection.bernoulli_head_V(i, j, k, grid, scheme, U.u, U.v)
    end
end

#####
##### Switching
#####

"""
    set_momentum_fix!(name)

`name` ∈ ("stock", "symdiv", "ke", "hybrid", "flux"). "ke" includes symdiv. "flux" uses no override (the model is
built with the flux form). Call before building the model, and run it through `Base.invokelatest`.
"""
function set_momentum_fix!(name)
    name ∈ ("stock", "symdiv", "ke", "hybrid", "flux") || error("unknown fix $name")
    set_symmetric_divergence!(name ∈ ("symdiv", "ke"))
    @eval ke_fix_on() = $(name == "ke")
    @eval hybrid_on() = $(name == "hybrid")
    return nothing
end

"""
    report_fix_faces(grid)

Count the wet u- and v-faces where "ke" and "hybrid" act (the cell thickness changes nearby), level by level summed.
A quick check, on the CPU, that the fixes reach the grid at all and touch only a small part of it.
"""
function report_fix_faces(grid)
    grid isa ImmersedBoundaryGrid || (println("not an immersed-boundary grid: the fixes do nothing here"); return nothing)
    cpu = Oceananigans.Architectures.on_architecture(Oceananigans.Architectures.CPU(), grid)
    Nx, Ny, Nz = size(cpu)
    C = Center(); F = Face()
    nu = au = nv = av = 0
    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        if !Oceananigans.Grids.peripheral_node(i, j, k, cpu, F, C, C)
            nu += 1; au += thickness_changes_x(i, j, k, cpu)
        end
        if acts_in_y(cpu) && !Oceananigans.Grids.peripheral_node(i, j, k, cpu, C, F, C)
            nv += 1; av += thickness_changes_y(i, j, k, cpu)
        end
    end
    println("fix acts at ", au, " of ", nu, " wet u-faces (", round(100au / max(nu, 1), digits = 1), " %)",
            acts_in_y(cpu) ? string(" and ", av, " of ", nv, " wet v-faces (", round(100av / max(nv, 1), digits = 1), " %)") : "")
    return nothing
end
