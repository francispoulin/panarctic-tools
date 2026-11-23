# Rotation cap on CATKE's stable length (rotation_cap.jl, 23 Nov), as an override inside Oceananigans:
#     ℓ = min(Cˢ d, √e/N, Cᶠ √e/|f|)     (stable_length_scaleᶜᶜᶠ / ᶜᶜᶜ, catke_mixing_length.jl l. 76, 92)
# Cᶠ and f come from a Ref, set per run with set_rotation_cap!(Cᶠ, f); Cᶠ = 0 or f = 0 is stock CATKE.
# Since ℓᶠ / ℓᴺ = Cᶠ N / |f|, the cap acts where N < |f| / Cᶠ: in and just below the mixed layer for Cᶠ ~ 0.1. In a
# neutral layer it is a Blackadar-type limit: with √e ≈ 2 u★, ℓ ≤ 2 Cᶠ u★/|f| (≈ 0.2 u★/|f| for Cᶠ = 0.1).
# CPU columns only (the Ref is read inside the kernels). Loaded once per process: redefining these methods again would
# invalidate compiled code (workflow.md).

let M = Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities
    isdefined(M, :rotation_cap) || Core.eval(M, quote
        const rotation_cap = Ref((Cᶠ = 0.0, f = 0.0))

        @inline function rotation_mixing_length(w★, FT)
            Cᶠ, f = rotation_cap[]
            off = (Cᶠ == 0) | (f == 0)
            return ifelse(off, FT(Inf), FT(Cᶠ * w★ / abs(f)))
        end

        @inline function stable_length_scaleᶜᶜᶠ(i, j, k, grid, closure, e, velocities, tracers, buoyancy)
            Cˢ = closure.mixing_length.Cˢ
            Cᵇ = closure.mixing_length.Cᵇ
            d_up   = Cˢ * depthᶜᶜᶠ(i, j, k, grid)
            d_down = Cᵇ * height_above_bottomᶜᶜᶠ(i, j, k, grid)
            d = min(d_up, d_down)
            ℓᴺ = stratification_mixing_lengthᶜᶜᶠ(i, j, k, grid, closure, e, tracers, buoyancy)
            w★ = ℑzᵃᵃᶠ(i, j, k, grid, turbulent_velocityᶜᶜᶜ, closure, e)
            ℓᶠ = rotation_mixing_length(w★, eltype(grid))
            ℓ = min(d, ℓᴺ, ℓᶠ)
            ℓ = ifelse(isnan(ℓ), d, ℓ)
            return ℓ
        end

        @inline function stable_length_scaleᶜᶜᶜ(i, j, k, grid, closure, e, velocities, tracers, buoyancy)
            Cˢ = closure.mixing_length.Cˢ
            Cᵇ = closure.mixing_length.Cᵇ
            d_up   = Cˢ * depthᶜᶜᶜ(i, j, k, grid)
            d_down = Cᵇ * height_above_bottomᶜᶜᶜ(i, j, k, grid)
            d = min(d_up, d_down)
            ℓᴺ = stratification_mixing_lengthᶜᶜᶜ(i, j, k, grid, closure, e, tracers, buoyancy)
            w★ = turbulent_velocityᶜᶜᶜ(i, j, k, grid, closure, e)
            ℓᶠ = rotation_mixing_length(w★, eltype(grid))
            ℓ = min(d, ℓᴺ, ℓᶠ)
            ℓ = ifelse(isnan(ℓ), d, ℓ)
            return ℓ
        end
    end)
end

set_rotation_cap!(Cᶠ, f) =
    (Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities.rotation_cap[] = (Cᶠ = Float64(Cᶠ), f = Float64(f)))
