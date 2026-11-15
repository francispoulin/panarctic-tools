# Thin partial bottom cells next to thick ones: bottom jet and extra mixing with WENOVectorInvariant.
# Self-contained example for the Oceananigans issue (plain Oceananigans, no overrides). CPU, about a minute.
#
# 2-D channel, periodic in x: a 2837 m basin with a plateau at 1863 m. With 20 exponential levels the plateau's bottom
# cell is 83 of 415 m wet. A current with the same transport everywhere crosses the plateau for 10 days. Nothing
# forces the flow, so the rise of the background potential energy (the PE of the same water sorted by density) is
# numerical mixing.

using Oceananigans, Oceananigans.Units, Printf
using Oceananigans.Grids: inactive_node, peripheral_node, znode
using Oceananigans.Operators: volume

const Nx, Nz, Δx = 32, 20, 15kilometers
const basin, plateau = 2837, 1863
on_plateau(x) = 12Δx ≤ x < 20Δx
bottom(x) = on_plateau(x) ? -plateau : -basin
uᵢ(x, z) = 0.1 * basin / (on_plateau(x) ? plateau : basin)     # same transport over basin and plateau
bᵢ(x, z) = 0.01 * exp(z / 1000)

function kinetic_energy(model)
    grid, u = model.grid, model.velocities.u
    return sum(u[i, 1, k]^2 / 2 * volume(i, 1, k, grid, Face(), Center(), Center())
               for i in 1:Nx, k in 1:Nz if !peripheral_node(i, 1, k, grid, Face(), Center(), Center()))
end

# Background potential energy: the water sorted by density into the model's own cells, heaviest lowest
function background_potential_energy(model)
    grid, b = model.grid, model.tracers.b
    C = Center()
    wet = [(i, k) for i in 1:Nx, k in 1:Nz if !inactive_node(i, 1, k, grid, C, C, C)]
    bs = [b[i, 1, k] for (i, k) in wet]
    zs = [znode(i, 1, k, grid, C, C, C) for (i, k) in wet]
    Vs = [volume(i, 1, k, grid, C, C, C) for (i, k) in wet]
    cells, parcels = sortperm(zs), sortperm(bs)
    BPE, c, room = 0.0, 1, Vs[cells[1]]
    for n in parcels
        left = Vs[n]
        while left > 1e-12 * Vs[n]
            take = min(left, room)
            BPE -= bs[n] * zs[cells[c]] * take
            left -= take; room -= take
            if room ≤ 1e-12 * Vs[cells[c]] && c < length(cells)
                c += 1; room = Vs[cells[c]]
            end
        end
    end
    return BPE
end

# Fastest flow in the lowest five wet u-faces of any column
function bottom_jet(model)
    grid, u = model.grid, model.velocities.u
    jet = 0.0
    for i in 1:Nx
        wet = [k for k in 1:Nz if !peripheral_node(i, 1, k, grid, Face(), Center(), Center())]
        for k in wet[1:min(5, end)]
            jet = max(jet, abs(u[i, 1, k]))
        end
    end
    return jet
end

function run_case(Bottom, momentum_advection; days = 10)
    z = ExponentialDiscretization(Nz, -5000, 0; scale = 1250, mutable = true)
    underlying = RectilinearGrid(size = (Nx, Nz), x = (0, Nx * Δx), z, halo = (7, 7),
                                 topology = (Periodic, Flat, Bounded))
    grid = ImmersedBoundaryGrid(underlying, Bottom(bottom))
    model = HydrostaticFreeSurfaceModel(grid; momentum_advection, tracer_advection = WENO(order = 5),
                                        buoyancy = BuoyancyTracer(), tracers = :b,
                                        free_surface = SplitExplicitFreeSurface(grid; substeps = 30))
    set!(model; u = uᵢ, b = bᵢ)
    KE₀, BPE₀ = kinetic_energy(model), background_potential_energy(model)
    jet = 0.0
    for day in 1:days
        for _ in 1:144
            time_step!(model, 10minutes)
        end
        jet = max(jet, bottom_jet(model))
    end
    return (background_potential_energy(model) - BPE₀) / KE₀, jet
end

println("mixing = rise of background PE / initial KE after 10 days; jet = fastest flow near the bottom")
@printf("%-28s %-12s %-14s %-8s %s\n", "momentum advection", "whole cells", "partial cells", "ratio", "jet, partial cells")
for (name, scheme) in (("WENOVectorInvariant()", WENOVectorInvariant()), ("WENO(order = 5), flux form", WENO(order = 5)))
    mw, _  = run_case(GridFittedBottom, scheme)
    mp, jp = run_case(PartialCellBottom, scheme)
    @printf("%-28s %-12.3f %-14.3f %-8.2f %.2f m/s\n", name, mw, mp, mp / mw, jp)
end
