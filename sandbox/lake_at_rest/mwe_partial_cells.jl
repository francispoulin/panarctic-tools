# Minimal example: a partial-cell bottom mixes more than a whole-cell bottom.
#
# 2-D (x, z), no rotation, periodic in x. A 2837 m basin with a 1863 m plateau (two 974 m cliffs), stratification
# b(z) = Δb exp(z/d), and a depth-uniform current with the same transport everywhere (0.1 m/s over the basin). No
# forcing, so the total energy E = KE + PE + E_η can only fall, except through mixing, which raises the background
# potential energy BPE (the PE of the same water re-sorted by buoyancy). Numerical mixing is the only mixing here.
#
# The test: after `days`, ΔBPE (divided by the initial KE) with partial cells should be no larger than with whole
# cells. In Oceananigans 0.117.2 it is larger, and the fastest flow sits in the partial bottom cell at the downstream
# cliff.
#
#   julia --project=<PanArcticShenanigans>/nested mwe_partial_cells.jl [days] [Nz] [zstar|static] [n_bottom] [plateau depth] [vector|centred|flux]
#
# The last argument picks the momentum advection: WENOVectorInvariant() ("vector", default), VectorInvariant()
# ("centred") or the flux form WENO(order = 5) ("flux").
#
# Consistency check at t = 0 (vector-invariant schemes): u starts uniform in depth in every column, so the vertical
# advection w ∂z u is exactly zero; the script prints the largest discrete value over all wet u-faces. This is the
# instant test for a fix.
#
# max|u| is over the lowest n_bottom wet u-faces of every column (default 5), so near-bottom flow is compared, not the
# surface. The plateau depth (default 1863 m) sets the wet fraction of the plateau's bottom cell.

using Oceananigans, Oceananigans.Units, Printf
using Oceananigans.Operators: volume, Azᶜᶜᶜ
using Oceananigans.Grids: peripheral_node, inactive_node, znode
using Oceananigans.Architectures: convert_to_device
using Oceananigans.Operators: Operators, V⁻¹ᶠᶜᶜ, δzᵃᵃᶜ, δxᶠᶜᶜ, Δx⁻¹ᶠᶜᶜ
using Oceananigans.Advection: Advection

days     = length(ARGS) ≥ 1 ? parse(Float64, ARGS[1]) : 10.0
Nz       = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 20
vertical = length(ARGS) ≥ 3 ? ARGS[3] : "zstar"
n_bottom = length(ARGS) ≥ 4 ? parse(Int, ARGS[4]) : 5      # speeds compared over the lowest n wet faces of each column
plateau  = length(ARGS) ≥ 5 ? parse(Float64, ARGS[5]) : 1863.0
momentum = length(ARGS) ≥ 6 ? ARGS[6] : "vector"
momentum ∈ ("vector", "centred", "flux") || error("momentum advection must be vector, centred or flux, got $momentum")
momentum_scheme = Dict("vector"  => WENOVectorInvariant(),   # WENO, upwinded divergence + flux-form vertical term
                       "centred" => VectorInvariant(),       # centred, energy-conserving w ∂z u (MITgcm's form)
                       "flux"    => WENO(order = 5))[momentum]

Nx, Δx = 32, 15kilometers
Lx     = Nx * Δx
depth  = 5000
basin = 2837
on_plateau(x) = 12Δx ≤ x < 20Δx
bottom(x) = on_plateau(x) ? -plateau : -basin

U, Δb, d = 0.1, 0.01, 1000
bᵢ(x, z) = Δb * exp(z / d)
uᵢ(x, z) = U * basin / (on_plateau(x) ? plateau : basin)

function cliff_model(bottom_kind; scheme = momentum_scheme)
    z = ExponentialDiscretization(Nz, -depth, 0; scale = depth / 4, mutable = vertical == "zstar")
    underlying = RectilinearGrid(CPU(); size = (Nx, Nz), x = (0, Lx), z, halo = (7, 7),
                                 topology = (Periodic, Flat, Bounded))
    Bottom = bottom_kind == :partial ? PartialCellBottom : GridFittedBottom
    grid = ImmersedBoundaryGrid(underlying, Bottom(bottom))
    model = HydrostaticFreeSurfaceModel(grid;
                                        momentum_advection = scheme,
                                        tracer_advection = WENO(order = 5),
                                        buoyancy = BuoyancyTracer(), tracers = :b,
                                        free_surface = SplitExplicitFreeSurface(grid; substeps = 30))
    set!(model; u = uᵢ, b = bᵢ)
    return model
end

# −Σ b z V over the wet cells, and its minimum over rearrangements of the water (heaviest parcels into the lowest
# cells, splitting parcels where cell volumes differ). z = znode, the same height the model's energetics use.
function potential_energies(model)
    grid, b = model.grid, model.tracers.b
    bs, zs, Vs = Float64[], Float64[], Float64[]
    for i in 1:Nx, k in 1:Nz
        inactive_node(i, 1, k, grid, Center(), Center(), Center()) && continue
        push!(bs, b[i, 1, k]); push!(zs, znode(i, 1, k, grid, Center(), Center(), Center()))
        push!(Vs, volume(i, 1, k, grid, Center(), Center(), Center()))
    end
    PE = -sum(bs .* zs .* Vs)
    cells, parcels = sortperm(zs), sortperm(bs)
    BPE, c, room = 0.0, 1, Vs[cells[1]]
    for n in parcels
        left = Vs[n]
        while left > 1e-12 * Vs[n]
            take = min(left, room)
            BPE -= bs[n] * zs[cells[c]] * take
            left -= take; room -= take
            if room ≤ 1e-12 * Vs[cells[c]]
                c == length(cells) && break
                c += 1; room = Vs[cells[c]]
            end
        end
    end
    return PE, BPE
end

function diagnostics(model)
    grid, u, η = model.grid, model.velocities.u, model.free_surface.displacement
    KE, umax, loc = 0.0, 0.0, (0, 0)
    for i in 1:Nx
        counted = 0
        for k in 1:Nz
            peripheral_node(i, 1, k, grid, Face(), Center(), Center()) && continue
            KE += u[i, 1, k]^2 / 2 * volume(i, 1, k, grid, Face(), Center(), Center())
            counted += 1
            counted ≤ n_bottom && abs(u[i, 1, k]) > umax && ((umax, loc) = (abs(u[i, 1, k]), (i, k)))
        end
    end
    Eη = sum(Oceananigans.defaults.gravitational_acceleration * η[i, 1, Nz+1]^2 / 2 * Azᶜᶜᶜ(i, 1, Nz+1, grid) for i in 1:Nx)
    PE, BPE = potential_energies(model)
    return (; KE, PE, BPE, Eη, umax, loc)
end

println("Flow over a cliff, $(Nz) levels, $(vertical), momentum advection: $(summary(momentum_scheme)). Energies divided by the initial KE.")
println("ΔBPE = numerical mixing. Cliffs at faces i = 13 and 21 (plateau cells 13–20), plateau at $(plateau) m.")
println("Bottom max|u| = largest speed over the lowest $(n_bottom) wet u-faces of each column.\n")

# Where the bottom is, level by level, for each kind of bottom: the wet thickness of each level in a basin column
# (i = 5) and a plateau column (i = 16). 0 = immersed.
let
    grids = Dict(kind => cliff_model(kind).grid for kind in (:whole, :partial))
    ug = grids[:whole].underlying_grid
    println("level   depth range (m)        basin column: whole  partial    plateau column: whole  partial")
    for k in 1:8
        zf = znode(1, 1, k, ug, Center(), Center(), Face()); zt = znode(1, 1, k + 1, ug, Center(), Center(), Face())
        wet(kind, i) = inactive_node(i, 1, k, grids[kind], Center(), Center(), Center()) ? 0.0 :
                       Oceananigans.Operators.Δzᶜᶜᶜ(i, 1, k, grids[kind])
        @printf("%3d   %7.0f to %7.0f      %16.0f %8.0f    %19.0f %8.0f%s\n", k, zf, zt,
                wet(:whole, 5), wet(:partial, 5), wet(:whole, 16), wet(:partial, 16), "")
    end
    println()
end
# t = 0: largest |vertical advection| over the wet u-faces; zero for a consistent scheme (u is uniform in depth)
if momentum != "flux"
    println("t = 0 check: u is uniform in depth, so the vertical advection w ∂z u must vanish. Largest |value| (m/s per day):")
    for bottom_kind in (:whole, :partial)
        model = cliff_model(bottom_kind)
        grid, Uk = convert_to_device(CPU(), (model.grid, model.velocities))
        worst, at = 0.0, (0, 0)
        for i in 1:Nx, k in 1:Nz
            peripheral_node(i, 1, k, model.grid, Face(), Center(), Center()) && continue
            a = abs(Advection.vertical_advection_U(i, 1, k, grid, model.advection.momentum, Uk))
            a > worst && ((worst, at) = (a, (i, k)))
        end
        @printf("  %-7s %10.4f at (i, k) = (%2d, %2d)\n", bottom_kind, 86400 * worst, at...)
    end
    println()
end

# Momentum budget at one u-face: each term's contribution to ∂u/∂t (m/s per day), as the model computes it, and the
# change of u over the next step. Vector-invariant form only (2-D, no rotation, so the vorticity term is ~0):
#   KE grad  = −∂x(u²/2)                       (bernoulli_head_U)
#   div      = −u ∇·U / V, upwinded            (upwinded_divergence_flux_U, part of vertical_advection_U)
#   ∂z(wu)   = −δz(W u) / V                    (vector_invariant_flux_Wu, the other part of vertical_advection_U)
#   p′ grad  = −∂x p′ (baroclinic, with the z-star chain rule)
#   g ∂x η   = −g δx η / Δx (barotropic; applied in the substeps, so not in "sum" but in "actual")
function momentum_budget(model, i, k)
    scheme = model.advection.momentum
    grid, U, p = convert_to_device(CPU(), (model.grid, model.velocities, model.pressure.pHY′))
    η = model.free_surface.displacement
    g = Oceananigans.defaults.gravitational_acceleration
    KE  = - Advection.bernoulli_head_U(i, 1, k, grid, scheme, U.u, U.v)
    hor = - Advection.horizontal_advection_U(i, 1, k, grid, scheme, U.u, U.v)
    div = - V⁻¹ᶠᶜᶜ(i, 1, k, grid) * Advection.upwinded_divergence_flux_Uᶠᶜᶜ(i, 1, k, grid, scheme, U.u, U.v)
    wu  = - V⁻¹ᶠᶜᶜ(i, 1, k, grid) * δzᵃᵃᶜ(i, 1, k, grid, Advection.vector_invariant_flux_Wu, scheme, U.w, U.u)
    pb  = - Operators.∂xᶠᶜᶜ(i, 1, k, grid, p)
    ps  = - g * (η[i, 1, Nz+1] - η[i-1, 1, Nz+1]) * Δx⁻¹ᶠᶜᶜ(i, 1, Nz+1, model.grid)
    return (; KE, hor, div, wu, pb, ps)
end

if momentum == "vector"
    hours = (0, 1, 2, 3, 6, 12, 24)
    for bottom_kind in (:whole, :partial)
        model = cliff_model(bottom_kind)
        i = 21                                    # downstream cliff
        k = findfirst(k -> !peripheral_node(i, 1, k, model.grid, Face(), Center(), Center()), 1:Nz)
        println("Momentum budget, $bottom_kind cells, face (i, k) = ($i, $k) (lowest wet face at the downstream cliff), m/s per day:")
        @printf("  %5s %8s %10s %10s %10s %10s %10s %10s %10s %10s\n", "hour", "u",
                "KE grad", "div", "∂z(wu)", "vorticity", "p′ grad", "sum", "actual", "g∂xη")
        nstep = 0
        for h in hours
            while nstep < 6h
                time_step!(model, 10minutes); nstep += 1
            end
            T = momentum_budget(model, i, k)
            u₀ = model.velocities.u[i, 1, k]
            time_step!(model, 10minutes); nstep += 1
            actual = (model.velocities.u[i, 1, k] - u₀) / 10minutes
            total = T.KE + T.div + T.wu + T.hor + T.pb
            @printf("  %5d %8.3f %+10.4f %+10.4f %+10.4f %+10.4f %+10.4f %+10.4f %+10.4f %+10.4f\n", h, u₀,
                    86400 .* (T.KE, T.div, T.wu, T.hor, T.pb, total, actual, T.ps)...)
        end
        println()
    end
end

# Same state, three discretisations of the advection term (u·∇)u. After one day of the chosen scheme (the jet is
# there), evaluate each form's contribution to ∂u/∂t (m/s per day) on that one state, at the faces around the
# downstream cliff, in the lowest two wet levels there. The flux form δx(Uu) + δz(Wu) moves momentum with the same
# transports and volumes as the water, so it is the reference; a vector-invariant piece that differs from it at the
# thin face is the inconsistency.
#   KE: −∂x(u²/2)    vert: the vertical term (w ∂z u)    flux x / flux z: −V⁻¹ δx(Uu), −V⁻¹ δz(Wu)
let
    for bottom_kind in (:whole, :partial)
        # Schemes as a model builds them for this grid (a bare WENOVectorInvariant() lacks the grid's settings)
        schemes = map(s -> cliff_model(bottom_kind; scheme = s).advection.momentum,
                      (; weno = WENOVectorInvariant(), centred = VectorInvariant(), flux = WENO(order = 5)))
        model = cliff_model(bottom_kind)
        for _ in 1:144
            time_step!(model, 10minutes)
        end
        grid, Uk = convert_to_device(CPU(), (model.grid, model.velocities))
        kb = findfirst(k -> !peripheral_node(21, 1, k, model.grid, Face(), Center(), Center()), 1:Nz)
        println("Advection term on the same state (after 1 day of $(momentum)), $bottom_kind cells, m/s per day:")
        @printf("  %8s %7s | %9s %9s %9s | %9s %9s %9s | %9s %9s %9s\n", "(i, k)", "u",
                "WENO-VI", "KE", "vert", "centred", "KE", "vert", "flux", "flux x", "flux z")
        for k in (kb, kb + 1), i in 19:23
            peripheral_node(i, 1, k, model.grid, Face(), Center(), Center()) && continue
            V⁻¹ = V⁻¹ᶠᶜᶜ(i, 1, k, grid)
            row = Float64[]
            for s in (schemes.weno, schemes.centred)
                KE   = - Advection.bernoulli_head_U(i, 1, k, grid, s, Uk.u, Uk.v)
                vert = - Advection.vertical_advection_U(i, 1, k, grid, s, Uk)
                hor  = - Advection.horizontal_advection_U(i, 1, k, grid, s, Uk.u, Uk.v)
                append!(row, (KE + vert + hor, KE, vert))
            end
            ftot = - Advection.div_𝐯u(i, 1, k, grid, schemes.flux, Uk, Uk.u)    # exactly what the flux-form model uses
            fz   = - V⁻¹ * δzᵃᵃᶜ(i, 1, k, grid, Advection._advective_momentum_flux_Wu, schemes.flux, Uk.w, Uk.u)
            append!(row, (ftot, ftot - fz, fz))
            @printf("  %8s %7.3f | %+9.4f %+9.4f %+9.4f | %+9.4f %+9.4f %+9.4f | %+9.4f %+9.4f %+9.4f\n",
                    "($i, $k)", model.velocities.u[i, 1, k], (86400 .* row)...)
        end
        println()
    end
end

results = Dict{Symbol, Any}()
for bottom_kind in (:whole, :partial)
    model = cliff_model(bottom_kind)
    D₀ = diagnostics(model)
    @printf("%-7s  %4s %9s %9s %9s %9s %16s\n", bottom_kind, "day", "ΔE", "ΔKE", "ΔBPE", "ΔAPE", "bottom max|u| (i, k)")
    local D = D₀
    for n in 1:ceil(Int, days)
        for _ in 1:144
            time_step!(model, 10minutes)
        end
        D = diagnostics(model)
        ΔE   = (D.KE + D.PE + D.Eη - D₀.KE - D₀.PE - D₀.Eη) / D₀.KE
        ΔBPE = (D.BPE - D₀.BPE) / D₀.KE
        ΔAPE = (D.PE - D.BPE - D₀.PE + D₀.BPE) / D₀.KE
        @printf("         %4d %+9.3f %+9.3f %+9.3f %+9.3f %12.3f (%2d,%2d)\n", n, ΔE, (D.KE - D₀.KE) / D₀.KE,
                ΔBPE, ΔAPE, D.umax, D.loc...)
    end
    results[bottom_kind] = (D.BPE - D₀.BPE) / D₀.KE
    println()
end

ratio = results[:partial] / results[:whole]
@printf("ΔBPE partial / whole = %.2f  →  %s\n", ratio, ratio ≤ 1.05 ? "PASS" : "FAIL (partial cells mix more)")
