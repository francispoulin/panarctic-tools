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
#   julia --project=<PanArcticShenanigans>/nested mwe_partial_cells.jl [days] [Nz] [zstar|static] [n_bottom] [plateau depth]
#
# max|u| is over the lowest n_bottom wet u-faces of every column (default 5), so near-bottom flow is compared, not the
# surface. The plateau depth (default 1863 m) sets the wet fraction of the plateau's bottom cell.

using Oceananigans, Oceananigans.Units, Printf
using Oceananigans.Operators: volume, Azᶜᶜᶜ
using Oceananigans.Grids: peripheral_node, inactive_node, znode

days     = length(ARGS) ≥ 1 ? parse(Float64, ARGS[1]) : 10.0
Nz       = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 20
vertical = length(ARGS) ≥ 3 ? ARGS[3] : "zstar"
n_bottom = length(ARGS) ≥ 4 ? parse(Int, ARGS[4]) : 5      # speeds compared over the lowest n wet faces of each column
plateau  = length(ARGS) ≥ 5 ? parse(Float64, ARGS[5]) : 1863.0

Nx, Δx = 32, 15kilometers
Lx     = Nx * Δx
depth  = 5000
basin = 2837
on_plateau(x) = 12Δx ≤ x < 20Δx
bottom(x) = on_plateau(x) ? -plateau : -basin

U, Δb, d = 0.1, 0.01, 1000
bᵢ(x, z) = Δb * exp(z / d)
uᵢ(x, z) = U * basin / (on_plateau(x) ? plateau : basin)

function cliff_model(bottom_kind)
    z = ExponentialDiscretization(Nz, -depth, 0; scale = depth / 4, mutable = vertical == "zstar")
    underlying = RectilinearGrid(CPU(); size = (Nx, Nz), x = (0, Lx), z, halo = (7, 7),
                                 topology = (Periodic, Flat, Bounded))
    Bottom = bottom_kind == :partial ? PartialCellBottom : GridFittedBottom
    grid = ImmersedBoundaryGrid(underlying, Bottom(bottom))
    model = HydrostaticFreeSurfaceModel(grid;
                                        momentum_advection = WENOVectorInvariant(),
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

println("Flow over a cliff, $(Nz) levels, $(vertical). Energies divided by the initial KE.")
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
