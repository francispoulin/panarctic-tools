# Shared test cases for the thin-bottom-cell fixes, used by pr_checks.jl.
#
# The same physical problem as compare_schemes.jl (a depth-uniform current over a plateau whose bottom cell is thin
# with partial cells; buoyancy b = Δb exp(z/d); nothing forces the flow, so any rise of the background potential
# energy is numerical mixing), in three geometries:
#
#   :x    2-D (x, z) channel, plateau across x (as compare_schemes.jl)
#   :y    the same channel turned by 90°: (y, z), the current is v. Exercises the v-face code of the fixes.
#   :xy   3-D, doubly periodic, a flat-topped circular seamount: plateau edges at every angle to the grid, and
#         staircase corners where the thickness changes in x and y at once.
#
# Scheme names: stock, symdiv, ke, hybrid, flux (as in thickness_jump_fixes.jl), and the same with "-wet" appended
# (e.g. "hybrid-wet"), which compares wet cells only when looking for a thickness change, so walls (cliff edges,
# coasts) do not count and whole-cell grids are left as stock.

include(joinpath(@__DIR__, "thickness_jump_fixes.jl"))

using Oceananigans.Units, Printf
using Oceananigans.Operators: volume, Azᶜᶜᶜ, Δrᶜᶜᶜ
using Oceananigans.Grids: peripheral_node, inactive_node, znode, znodes
using Oceananigans.Architectures: on_architecture

const Δh = 15kilometers
const depth, basin, default_plateau = 5000, 2837, 1863.0
const Uc, Δb, dscale = 0.1, 0.01, 1000
const n_bottom = 5                       # "near-bottom" = the lowest 5 wet faces of each column

Base.@kwdef struct ThinCellCase
    geometry :: Symbol  = :x             # :x, :y or :xy
    Nz       :: Int     = 20
    plateau  :: Float64 = default_plateau
    vertical :: Symbol  = :zstar         # :zstar or :static
    latitude :: Float64 = 0.0            # 0: no rotation
    arch               = CPU()
end

function with(c::ThinCellCase; kw...)
    d = Dict{Symbol, Any}(f => getfield(c, f) for f in fieldnames(ThinCellCase))
    for (k, v) in kw
        d[k] = v
    end
    return ThinCellCase(; d...)
end

describe(c::ThinCellCase) = @sprintf("%s, %d levels, %s, %s, plateau %.1f m", c.geometry == :xy ? "3-D seamount" :
                                     "2-D channel in $(c.geometry)", c.Nz, c.vertical == :zstar ? "z-star" : "static",
                                     c.latitude == 0 ? "no rotation" : "f at $(c.latitude)°N", c.plateau)

horizontal_size(c) = c.geometry == :xy ? 40 : 32
on_plateau_1d(s) = 12Δh ≤ s < 20Δh
on_seamount(x, y) = (x - 20Δh)^2 + (y - 20Δh)^2 < (8Δh)^2

function underlying_grid(c)
    N = horizontal_size(c)
    L = N * Δh
    z = ExponentialDiscretization(c.Nz, -depth, 0; scale = depth / 4, mutable = c.vertical == :zstar)
    if c.geometry == :x
        return RectilinearGrid(c.arch; size = (N, c.Nz), x = (0, L), z, halo = (7, 7),
                               topology = (Periodic, Flat, Bounded))
    elseif c.geometry == :y
        return RectilinearGrid(c.arch; size = (N, c.Nz), y = (0, L), z, halo = (7, 7),
                               topology = (Flat, Periodic, Bounded))
    else
        return RectilinearGrid(c.arch; size = (N, N, c.Nz), x = (0, L), y = (0, L), z, halo = (7, 7, 7),
                               topology = (Periodic, Periodic, Bounded))
    end
end

function bottom_function(c)
    P = c.plateau
    c.geometry == :x && return x -> on_plateau_1d(x) ? -P : -basin
    c.geometry == :y && return y -> on_plateau_1d(y) ? -P : -basin
    return (x, y) -> on_seamount(x, y) ? -P : -basin
end

# The same transport everywhere: H u = U × basin depth (divergence-free, as in compare_schemes.jl)
function initial_conditions(c)
    P = c.plateau
    if c.geometry == :x
        return (; u = (x, z) -> Uc * basin / (on_plateau_1d(x) ? P : basin), b = (x, z) -> Δb * exp(z / dscale))
    elseif c.geometry == :y
        return (; v = (y, z) -> Uc * basin / (on_plateau_1d(y) ? P : basin), b = (y, z) -> Δb * exp(z / dscale))
    else
        return (; u = (x, y, z) -> Uc * basin / (on_seamount(x, y) ? P : basin), b = (x, y, z) -> Δb * exp(z / dscale))
    end
end

function thin_cell_grid(c, bottom_kind; bottom = bottom_function(c))
    Bottom = bottom_kind == :partial ? PartialCellBottom : GridFittedBottom
    return ImmersedBoundaryGrid(underlying_grid(c), Bottom(bottom))
end

function thin_cell_model(c, bottom_kind, base_scheme; bottom = bottom_function(c))
    grid = thin_cell_grid(c, bottom_kind; bottom)
    momentum_advection = base_scheme == "flux" ? WENO(order = 5) : WENOVectorInvariant()
    coriolis = c.latitude == 0 ? nothing : FPlane(latitude = c.latitude)
    model = HydrostaticFreeSurfaceModel(grid; momentum_advection, coriolis,
                                        tracer_advection = WENO(order = 5),
                                        buoyancy = BuoyancyTracer(), tracers = :b,
                                        free_surface = SplitExplicitFreeSurface(grid; substeps = 30))
    set!(model; initial_conditions(c)...)
    return model
end

# Set the override switches for a scheme name; returns the base name ("hybrid-wet" → "hybrid")
function use_scheme!(name)
    base = replace(name, "-wet" => "")
    set_momentum_fix!(base)
    set_jump_walls!(!endswith(name, "-wet"))
    return base
end

# Wet fraction of the plateau's bottom cell with partial cells
function wet_fraction(c)
    grid = thin_cell_grid(with(c; arch = CPU()), :partial)
    i, j = c.geometry == :x ? (16, 1) : c.geometry == :y ? (1, 16) : (20, 20)
    kb = findfirst(k -> !inactive_node(i, j, k, grid, Center(), Center(), Center()), 1:c.Nz)
    return Δrᶜᶜᶜ(i, j, kb, grid) / Δrᶜᶜᶜ(i, j, kb, grid.underlying_grid)
end

# Plateau depth whose bottom cell is a fraction φ wet: the cell that holds the default plateau (1863 m), floor at
# φ of its thickness below its top face.
function plateau_for_fraction(φ, c)
    zf = znodes(underlying_grid(with(c; arch = CPU(), vertical = :static)), Face())   # same reference faces as z-star
    k = findlast(z -> z ≤ -default_plateau, zf)
    z_bottom, z_top = zf[k], zf[k+1]
    return -(z_top - φ * (z_top - z_bottom))
end

#####
##### Measurements (on CPU copies)
#####

cpu(x) = on_architecture(CPU(), x)

function potential_energies(grid, b)
    Nx, Ny, Nz = size(grid)
    C = Center()
    bs, zs, Vs = Float64[], Float64[], Float64[]
    for j in 1:Ny, i in 1:Nx, k in 1:Nz
        inactive_node(i, j, k, grid, C, C, C) && continue
        push!(bs, b[i, j, k]); push!(zs, znode(i, j, k, grid, C, C, C))
        push!(Vs, volume(i, j, k, grid, C, C, C))
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

# Kinetic, potential (and background) and free-surface energy; the near-bottom jet is the largest in-plane speed in
# the lowest n_bottom wet faces of any column. In a Flat direction the out-of-plane velocity sits at cell centres.
function measure(model, c)
    grid = cpu(model.grid)
    u, v = cpu(model.velocities.u), cpu(model.velocities.v)
    b, η = cpu(model.tracers.b), cpu(model.free_surface.displacement)
    Nx, Ny, Nz = size(grid)
    C, F = Center(), Face()
    ℓu = c.geometry == :y ? (C, C, C) : (F, C, C)
    ℓv = c.geometry == :x ? (C, C, C) : (C, F, C)
    jet_u, jet_v = c.geometry != :y, c.geometry != :x
    KE, jet = 0.0, 0.0
    for j in 1:Ny, i in 1:Nx
        nu = nv = 0
        for k in 1:Nz
            if !peripheral_node(i, j, k, grid, ℓu...)
                KE += u[i, j, k]^2 / 2 * volume(i, j, k, grid, ℓu...)
                nu += 1
                jet_u && nu ≤ n_bottom && (jet = max(jet, abs(u[i, j, k])))
            end
            if !peripheral_node(i, j, k, grid, ℓv...)
                KE += v[i, j, k]^2 / 2 * volume(i, j, k, grid, ℓv...)
                nv += 1
                jet_v && nv ≤ n_bottom && (jet = max(jet, abs(v[i, j, k])))
            end
        end
    end
    g = Oceananigans.defaults.gravitational_acceleration
    Eη = sum(g * η[i, j, Nz+1]^2 / 2 * Azᶜᶜᶜ(i, j, Nz+1, grid) for j in 1:Ny, i in 1:Nx)
    PE, BPE = potential_energies(grid, b)
    return (; KE, PE, BPE, Eη, jet)
end

#####
##### Runs
#####

function _run_case(c, bottom_kind, base, days)
    model = thin_cell_model(c, bottom_kind, base)
    E₀ = measure(model, c)
    E, jet = E₀, 0.0
    for day in 1:ceil(Int, days)
        for _ in 1:144
            time_step!(model, 10minutes)
        end
        E = measure(model, c)
        isfinite(E.KE) || break
        jet = max(jet, E.jet)
    end
    mixing = (E.BPE - E₀.BPE) / E₀.KE
    energy = (E.KE + E.PE + E.Eη - E₀.KE - E₀.PE - E₀.Eη) / E₀.KE
    return (; mixing, energy, jet)
end

"""
    run_case(c, bottom_kind, name; days = 10)

One run (`bottom_kind` :whole or :partial) with scheme `name`; returns mixing (ΔBPE / KE₀), energy (ΔE / KE₀) and
the largest daily near-bottom speed.
"""
function run_case(c, bottom_kind, name; days = 10)
    base = use_scheme!(name)
    return Base.invokelatest(_run_case, c, bottom_kind, base, days)
end

"""
    scheme_table(c, schemes; days = 10)

Whole and partial cells for each scheme, and the compare_schemes.jl verdict: partial mixes ≤ 1.2 × whole, the
partial-cell jet ≤ 1.1 × the flux form's, whole cells within 15 % of stock, and the energy falls.
"""
function scheme_table(c, schemes; days = 10)
    @printf("\n%s. Plateau bottom cell %.0f %% wet with partial cells. %g days.\n", describe(c), 100wet_fraction(c), days)
    results = Dict{Tuple{String, Symbol}, Any}()
    for name in schemes, bottom_kind in (:whole, :partial)
        t = @elapsed results[(name, bottom_kind)] = run_case(c, bottom_kind, name; days)
        @printf("  done: %-11s %-8s (%.0f s)\n", name, bottom_kind, t)
    end
    stock_whole  = haskey(results, ("stock", :whole)) ? results[("stock", :whole)].mixing : NaN
    flux_partial = haskey(results, ("flux", :partial)) ? results[("flux", :partial)].jet : NaN
    println("scheme       mixing (whole → partial)   ratio   jet, partial   whole vs stock   energy, partial    verdict")
    println("-----------  ------------------------  ------  -------------  ---------------  -----------------  -------")
    for name in schemes
        w, p = results[(name, :whole)], results[(name, :partial)]
        ratio = p.mixing / w.mixing
        change = isnan(stock_whole) ? NaN : 100 * (w.mixing / stock_whole - 1)
        whole_kept = name == "flux" || isnan(change) || abs(change) ≤ 15
        jet_ok = isnan(flux_partial) || p.jet ≤ 1.1 * flux_partial
        pass = ratio ≤ 1.2 && jet_ok && whole_kept && p.energy ≤ 0
        note = name == "flux" ? "reference" : name == "stock" ? (pass ? "PASS" : "FAIL (expected)") : pass ? "PASS" : "FAIL"
        @printf("%-11s  %8.3f → %-8.3f         %5.2f   %6.3f m/s     %+7.2f %%        %+7.3f %-8s  %s\n",
                name, w.mixing, p.mixing, ratio, p.jet, change, p.energy, p.energy > 0 ? "(rises)" : "(falls)", note)
    end
    return results
end
