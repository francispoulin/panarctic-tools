# One table: does a momentum-advection scheme stop thin bottom cells from mixing and making a bottom jet?
#
#   julia --project=<PanArcticShenanigans>/nested compare_schemes.jl [days] [latitude] [schemes...]
#
#   days      default 10
#   latitude  default 0 (no rotation); e.g. 72 for f at 72°N
#   schemes   default: stock symdiv ke hybrid flux
#
# The test (same as mwe_partial_cells.jl, without its detailed diagnostics): a 2-D channel, 2837 m deep, with a
# 1863 m plateau whose bottom cell is only 20 % wet with partial cells. A depth-uniform current crosses the plateau
# for `days`. Nothing forces the flow, so any rise of the background potential energy is numerical mixing.
# Each scheme runs twice: whole cells and partial cells. CPU; a few minutes per scheme.
#
# Schemes:
#   stock   WENOVectorInvariant(), as the pan-Arctic model uses now
#   symdiv  stock + the u-cell's own divergence in the vertical term          (symmetric_divergence_fix.jl)
#   ke      symdiv + the flux form's u∂xu where the cell thickness changes    (thickness_jump_fixes.jl)
#   hybrid  stock, but the flux form where the cell thickness changes         (thickness_jump_fixes.jl)
#   flux    WENO(order = 5) flux form everywhere: the reference, as in MITgcm's default

include(joinpath(@__DIR__, "symmetric_divergence_fix.jl"))
include(joinpath(@__DIR__, "thickness_jump_fixes.jl"))

using Oceananigans.Units, Printf
using Oceananigans.Operators: volume, Azᶜᶜᶜ, Δrᶜᶜᶜ
using Oceananigans.Grids: peripheral_node, inactive_node, znode
using Oceananigans.Architectures: convert_to_device

days     = length(ARGS) ≥ 1 ? parse(Float64, ARGS[1]) : 10.0
latitude = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 0.0
schemes  = length(ARGS) ≥ 3 ? ARGS[3:end] : ["stock", "symdiv", "ke", "hybrid", "flux"]

const Nx, Nz, Δx = 32, 20, 15kilometers
const depth, basin, plateau = 5000, 2837, 1863
const U, Δb, d = 0.1, 0.01, 1000
const n_bottom = 5                               # "near-bottom" = the lowest 5 wet u-faces of each column

on_plateau(x) = 12Δx ≤ x < 20Δx
bottom(x) = on_plateau(x) ? -plateau : -basin
bᵢ(x, z) = Δb * exp(z / d)
uᵢ(x, z) = U * basin / (on_plateau(x) ? plateau : basin)

function cliff_model(bottom_kind, name)
    z = ExponentialDiscretization(Nz, -depth, 0; scale = depth / 4, mutable = true)
    underlying = RectilinearGrid(CPU(); size = (Nx, Nz), x = (0, Nx * Δx), z, halo = (7, 7),
                                 topology = (Periodic, Flat, Bounded))
    grid = ImmersedBoundaryGrid(underlying, bottom_kind == :partial ? PartialCellBottom(bottom) : GridFittedBottom(bottom))
    momentum_advection = name == "flux" ? WENO(order = 5) : WENOVectorInvariant()
    coriolis = latitude == 0 ? nothing : FPlane(; latitude)
    model = HydrostaticFreeSurfaceModel(grid; momentum_advection, coriolis,
                                        tracer_advection = WENO(order = 5),
                                        buoyancy = BuoyancyTracer(), tracers = :b,
                                        free_surface = SplitExplicitFreeSurface(grid; substeps = 30))
    set!(model; u = uᵢ, b = bᵢ)
    return model
end

# Potential energy −Σ b z V, and its minimum over rearrangements of the water (the background PE)
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

function energies(model)
    grid, u, v, η = model.grid, model.velocities.u, model.velocities.v, model.free_surface.displacement
    KE, jet = 0.0, 0.0
    for i in 1:Nx
        counted = 0
        for k in 1:Nz
            peripheral_node(i, 1, k, grid, Face(), Center(), Center()) && continue
            KE += (u[i, 1, k]^2 + v[i, 1, k]^2) / 2 * volume(i, 1, k, grid, Face(), Center(), Center())
            counted += 1
            counted ≤ n_bottom && (jet = max(jet, abs(u[i, 1, k])))
        end
    end
    g = Oceananigans.defaults.gravitational_acceleration
    Eη = sum(g * η[i, 1, Nz+1]^2 / 2 * Azᶜᶜᶜ(i, 1, Nz+1, grid) for i in 1:Nx)
    PE, BPE = potential_energies(model)
    return (; KE, PE, BPE, Eη, jet)
end

# Net push of the advection term alone on the domain's total x-momentum, Σ V (u·∇)u, as a percentage of Σ V |u| per
# day. A scheme that conserves momentum gives ~0 (its fluxes telescope); pressure and bottom drag are not included,
# so the real form drag of the cliffs does not enter.
function momentum_drift(model)
    grid, Uk = convert_to_device(CPU(), (model.grid, model.velocities))
    net, total = 0.0, 0.0
    for i in 1:Nx, k in 1:Nz
        peripheral_node(i, 1, k, model.grid, Face(), Center(), Center()) && continue
        V = volume(i, 1, k, model.grid, Face(), Center(), Center())
        net   += V * Advection.U_dot_∇u(i, 1, k, grid, model.advection.momentum, Uk)
        total += V * abs(Uk.u[i, 1, k])
    end
    return 100 * 86400 * net / total
end

# Run one case; the jet is the largest near-bottom speed, sampled once a day
function run_case_daily(bottom_kind, name)
    model = cliff_model(bottom_kind, name)
    E₀ = energies(model)
    jet, E = 0.0, E₀
    for day in 1:ceil(Int, days)
        for _ in 1:144
            time_step!(model, 10minutes)
        end
        E = energies(model)
        jet = max(jet, E.jet)
    end
    mixing = (E.BPE - E₀.BPE) / E₀.KE
    drift = momentum_drift(model)
    energy = (E.KE + E.PE + E.Eη - E₀.KE - E₀.PE - E₀.Eη) / E₀.KE
    return (; mixing, energy, jet, drift)
end

# Wet fraction of the plateau's bottom cell (partial cells)
φ = let grid = cliff_model(:partial, "stock").grid
    kb = findfirst(k -> !inactive_node(16, 1, k, grid, Center(), Center(), Center()), 1:Nz)
    Δrᶜᶜᶜ(16, 1, kb, grid) / Δrᶜᶜᶜ(16, 1, kb, grid.underlying_grid)
end

@printf("Thin bottom cells: %g days, %d levels, %s. Plateau bottom cell %.0f %% wet with partial cells.\n",
        days, Nz, latitude == 0 ? "no rotation" : "f at $(latitude)°N", 100φ)
println("Running: ", join(schemes, ", "), " (whole and partial cells each)\n")

results = Dict{Tuple{String, Symbol}, Any}()
for name in schemes, bottom_kind in (:whole, :partial)
    set_momentum_fix!(name)
    t = @elapsed results[(name, bottom_kind)] = Base.invokelatest(run_case_daily, bottom_kind, name)
    @printf("  done: %-7s %-8s (%.0f s)\n", name, bottom_kind, t)
end

# Verdict: partial cells must mix no more than 1.2 × whole cells, the near-bottom jet with partial cells must be no
# faster than 1.15 × the whole-cell jet, and (for a fix, not the flux reference) whole cells must stay within 15 %
# of the stock scheme's whole-cell mixing (the fix should not change the model where nothing is wrong).
stock_whole  = haskey(results, ("stock", :whole)) ? results[("stock", :whole)].mixing : NaN
flux_partial = haskey(results, ("flux", :partial)) ? results[("flux", :partial)].jet : NaN
println("scheme   mixing (whole → partial)   ratio   jet, partial   whole vs stock   energy, partial   momentum push   verdict")
println("-------  ------------------------  ------  -------------  ---------------  ----------------  --------------  -------")
for name in schemes
    w, p = results[(name, :whole)], results[(name, :partial)]
    ratio = p.mixing / w.mixing
    change = isnan(stock_whole) ? NaN : 100 * (w.mixing / stock_whole - 1)
    whole_kept = name == "flux" || isnan(change) || abs(change) ≤ 15
    jet_ok = isnan(flux_partial) ? true : p.jet ≤ 1.1 * flux_partial
    pass = ratio ≤ 1.2 && jet_ok && whole_kept && p.energy ≤ 0
    note = name == "flux" ? "reference" : pass ? "PASS" : "FAIL"
    @printf("%-7s  %8.3f → %-8.3f         %5.2f   %6.2f m/s     %+6.0f %%          %+7.3f %-8s  %+7.3f %%/day   %s\n",
            name, w.mixing, p.mixing, ratio, p.jet, change, p.energy, p.energy > 0 ? "(rises)" : "(falls)",
            p.drift, note)
end
println("""

How to read it
  mixing          rise of the background potential energy over the run, in units of the initial kinetic energy.
                  Nothing forces the flow, so this is numerical mixing. Lower is better.
  ratio           partial-cell mixing ÷ whole-cell mixing. 1 = thin bottom cells add nothing; stock gives ~2.4.
                  Pass: 1.2 or less.
  jet, partial    fastest flow in the lowest 5 wet faces of any column, partial cells (largest daily value).
                  Pass: within 10 % of the flux form's.
  whole vs stock  change of the whole-cell mixing from the stock scheme. A fix should not change the model where
                  nothing was wrong. Pass: within 15 %.
  energy          change of total energy (kinetic + potential + free surface), partial cells. With no forcing it
                  must fall; a rise means the scheme feeds the flow. Pass: falls.
  momentum push   net push of the advection term alone on the channel's total momentum at the end of the run,
                  % of the momentum per day, partial cells. ~0 = momentum conserved by advection (flux form).
                  Shown, not part of the verdict: the vector-invariant form never conserves it exactly.
  "flux" is the reference (MITgcm's default form), not a candidate.""")
