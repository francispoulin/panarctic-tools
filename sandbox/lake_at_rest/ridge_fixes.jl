# Flow over a ridge with the experimental bottom-cell fixes (partial_cell_fixes.jl), and a mixing/conversion split of
# the potential-energy change.
#
#   julia --project=<PanArcticShenanigans>/nested ridge_fixes.jl [profile.csv] [days] [Nz]
#
# Same set-up as flow_over_ridge.jl (periodic channel, two 1000 m cliffs, z-star, GLORYS profile, linear equation of
# state, barotropic 0.1 m/s current, f = 0 and f(72°N)), run for these cases:
#
#   partial cells: none, F1, F2, F1+F2          whole cells: none, F1, F1+F2
#   (F2 alone changes nothing: see partial_cell_fixes.jl)
#
# "none" is the stock model and should repeat flow_over_ridge.jl (day 10, Nz = 20: ΔE/KE₀ ≈ +0.26 partial vs +0.12
# whole at f = 0; at f(72°N) partial gains about +0.09 from day 2 to 10, whole +0.015).
#
# Energies (per unit ρ₀, all divided by the initial kinetic energy KE₀):
#   PE  = −Σ b z V with z = znode (full-cell centres, as the model's own energetics use);
#   BPE = the same sum after re-sorting the water: parcels ordered by b, heaviest placed in the lowest cells (by znode),
#         filling cell volumes in order. BPE is the smallest PE any rearrangement can have, so it can only rise
#         through mixing. ΔBPE = mixing; ΔAPE = ΔPE − ΔBPE = the reversible exchange with kinetic energy.
#   A spurious source shows up in ΔAPE (and ΔE); excess numerical mixing shows up in ΔBPE.
#
# One Julia session runs every case (the fixes are switched, not reloaded). CPU only.

using Oceananigans, Oceananigans.Units, Printf
using Oceananigans.Operators: volume, Azᶜᶜᶜ
using Oceananigans.Grids: peripheral_node, inactive_node, znode, Center
using Oceananigans.BuoyancyFormulations: buoyancy_perturbationᶜᶜᶜ

include(joinpath(@__DIR__, "partial_cell_fixes.jl"))

profile_file = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(@__DIR__, "glorys_profile_0E_72N.csv")
run_days     = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 10.0
Nz_run       = length(ARGS) ≥ 3 ? parse(Int, ARGS[3]) : 20

rows = [parse.(Float64, split(line, ',')) for line in readlines(profile_file)[2:end] if !isempty(line)]
profile_depth = first.(rows)
profile_T     = getindex.(rows, 2)
profile_S     = getindex.(rows, 3)

function interpolate_profile(values, z)
    d = clamp(-z, first(profile_depth), last(profile_depth))
    k = clamp(searchsortedlast(profile_depth, d), 1, length(profile_depth) - 1)
    α = (d - profile_depth[k]) / (profile_depth[k + 1] - profile_depth[k])
    return (1 - α) * values[k] + α * values[k + 1]
end

Tᵢ(x, y, z) = interpolate_profile(profile_T, z)
Sᵢ(x, y, z) = interpolate_profile(profile_S, z)

Nx, Ny  = 32, 8
spacing = 15kilometers
Lx, Ly  = Nx * spacing, Ny * spacing
depth   = 5000
basin_depth, plateau_depth = 2837, 1863
on_plateau(x) = 12spacing ≤ x < 20spacing
bottom(x, y) = on_plateau(x) ? -plateau_depth : -basin_depth

U = 0.1
g = Oceananigans.defaults.gravitational_acceleration
uᵢ(x, y, z) = U * basin_depth / (on_plateau(x) ? plateau_depth : basin_depth)

function ridge_model(bottom_kind, Nz, latitude)
    z = ExponentialDiscretization(Nz, -depth, 0; scale = depth / 4, mutable = true)
    underlying_grid = RectilinearGrid(CPU(); size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z,
                                      halo = (7, 7, 7), topology = (Periodic, Bounded, Bounded))
    Bottom = bottom_kind == :partial ? PartialCellBottom : GridFittedBottom
    grid = ImmersedBoundaryGrid(underlying_grid, Bottom(bottom))

    coriolis = isnothing(latitude) ? nothing : FPlane(; latitude)
    f = isnothing(latitude) ? 0.0 : coriolis.f

    model = HydrostaticFreeSurfaceModel(grid;
                                        momentum_advection = WENOVectorInvariant(),
                                        tracer_advection = WENO(order = 5),
                                        buoyancy = SeawaterBuoyancy(equation_of_state = LinearEquationOfState(thermal_expansion = 2e-4, haline_contraction = 7.8e-4)),
                                        coriolis,
                                        free_surface = SplitExplicitFreeSurface(grid; substeps = 30),
                                        tracers = (:T, :S))

    set!(model; η = (x, y, z...) -> -f * U * (y - Ly / 2) / g)
    set!(model; u = uᵢ, T = Tᵢ, S = Sᵢ)
    return model
end

bottom_face(i, j, k, grid) = inactive_node(i - 1, j, k - 1, grid, Center(), Center(), Center()) |
                             inactive_node(i,     j, k - 1, grid, Center(), Center(), Center())

# Smallest −Σ b z V over rearrangements of the water among the model's cells: heaviest parcels (lowest b) fill the
# lowest cells first, splitting parcels where cell volumes differ.
function background_potential_energy(b, z, V)
    cells   = sortperm(z)
    parcels = sortperm(b)
    BPE  = 0.0
    c    = 1
    room = V[cells[c]]
    for n in parcels
        left = V[n]
        while left > 1e-12 * V[n]
            take = min(left, room)
            BPE -= b[n] * z[cells[c]] * take
            left -= take
            room -= take
            if room ≤ 1e-12 * V[cells[c]]
                c == length(cells) && break
                c += 1
                room = V[cells[c]]
            end
        end
    end
    return BPE
end

function energies(model)
    grid = model.grid
    Nz = size(grid, 3)
    u, v, η = model.velocities.u, model.velocities.v, model.free_surface.displacement
    kinetic = bottom_kinetic = potential = surface = 0.0
    fastest = (0.0, 0, 0, false)
    bs, zs, Vs = Float64[], Float64[], Float64[]
    for i in 1:Nx, j in 1:Ny, k in 1:Nz
        if !peripheral_node(i, j, k, grid, Face(), Center(), Center())
            e = u[i, j, k]^2 / 2 * volume(i, j, k, grid, Face(), Center(), Center())
            kinetic += e
            at_bottom = bottom_face(i, j, k, grid)
            at_bottom && (bottom_kinetic += e)
            abs(u[i, j, k]) > abs(fastest[1]) && (fastest = (u[i, j, k], i, k, at_bottom))
        end
        peripheral_node(i, j, k, grid, Center(), Face(), Center()) ||
            (kinetic += v[i, j, k]^2 / 2 * volume(i, j, k, grid, Center(), Face(), Center()))
        if !inactive_node(i, j, k, grid, Center(), Center(), Center())
            b = buoyancy_perturbationᶜᶜᶜ(i, j, k, grid, model.buoyancy.formulation, model.tracers)
            z = znode(i, j, k, grid, Center(), Center(), Center())
            V = volume(i, j, k, grid, Center(), Center(), Center())
            potential -= b * z * V
            push!(bs, b); push!(zs, z); push!(Vs, V)
        end
    end
    for i in 1:Nx, j in 1:Ny
        surface += g * η[i, j, Nz + 1]^2 / 2 * Azᶜᶜᶜ(i, j, Nz + 1, grid)
    end
    background = background_potential_energy(bs, zs, Vs)
    return (; kinetic, bottom_kinetic, potential, background, surface, fastest)
end

# Set RIDGE_NO_RUN = true before `include` to get the definitions only (probe_fixes.jl does).
if !(isdefined(Main, :RIDGE_NO_RUN) && Main.RIDGE_NO_RUN)

cases = [(:partial, (F1 = false, F2 = false)),
         (:partial, (F1 = true,  F2 = false)),
         (:partial, (F1 = false, F2 = true)),
         (:partial, (F1 = true,  F2 = true)),
         (:whole,   (F1 = false, F2 = false)),
         (:whole,   (F1 = true,  F2 = false)),
         (:whole,   (F1 = true,  F2 = true))]

results = []

println("Energy changes divided by the initial kinetic energy KE₀. ΔAPE = ΔPE − ΔBPE (conversion), ΔBPE = mixing.")
println("max|u| in m/s at (i, k), * = a bottom face. The cliffs are the faces i = 13 and i = 21. Nz = $Nz_run.\n")

for latitude in (nothing, 72), (bottom_kind, fixes) in cases
    set_fixes!(; fixes...)
    model = ridge_model(bottom_kind, Nz_run, latitude)
    E₀ = energies(model)
    KE₀ = E₀.kinetic
    total₀ = E₀.kinetic + E₀.potential + E₀.surface
    flabel = isnothing(latitude) ? "f = 0" : "f(72°N)"
    label = @sprintf("%-7s cells, %-5s, %s", bottom_kind, fix_label(), flabel)
    println(label)
    @printf("  %4s %9s %9s %9s %9s %9s %9s %9s %16s\n",
            "day", "ΔE", "ΔKE", "ΔKE_bot", "ΔPE", "ΔBPE", "ΔAPE", "ΔE_η", "max|u| (i, k)")
    Δt = 10minutes
    row = nothing
    for n in 1:ceil(Int, run_days)
        for _ in 1:round(Int, Oceananigans.Units.day / Δt)
            time_step!(model, Δt)
        end
        E = energies(model)
        u, i, k, at_bottom = E.fastest
        ΔE   = (E.kinetic + E.potential + E.surface - total₀) / KE₀
        ΔKE  = (E.kinetic - KE₀) / KE₀
        ΔPE  = (E.potential - E₀.potential) / KE₀
        ΔBPE = (E.background - E₀.background) / KE₀
        ΔAPE = ΔPE - ΔBPE
        Δη   = (E.surface - E₀.surface) / KE₀
        @printf("  %4d %+9.2e %+9.2e %+9.2e %+9.2e %+9.2e %+9.2e %+9.2e %7.3f (%2d,%2d)%s\n", n,
                ΔE, ΔKE, (E.bottom_kinetic - E₀.bottom_kinetic) / KE₀, ΔPE, ΔBPE, ΔAPE, Δη,
                abs(u), i, k, at_bottom ? "*" : "")
        row = (; label, ΔE, ΔKE, ΔPE, ΔBPE, ΔAPE, umax = abs(u), at_bottom)
    end
    FIX_F1[] && println("  F1 override calls: ", F1_CALLS[], F1_CALLS[] == 0 ? "  <-- F1 NOT IN USE" : "")
    println()
    push!(results, row)
end

set_fixes!()

println("Summary at day $(ceil(Int, run_days)), Nz = $Nz_run (energies / KE₀):")
@printf("  %-30s %9s %9s %9s %9s %9s %8s\n", "case", "ΔE", "ΔKE", "ΔPE", "ΔBPE", "ΔAPE", "max|u|")
for r in results
    @printf("  %-30s %+9.3f %+9.3f %+9.3f %+9.3f %+9.3f %7.3f%s\n",
            r.label, r.ΔE, r.ΔKE, r.ΔPE, r.ΔBPE, r.ΔAPE, r.umax, r.at_bottom ? "*" : "")
end

end
