# Stage 3: the pan-Arctic ocean and sea ice nested in GLORYS (NumericalEarth PR #676, ss/nested-ocean-model), forced by
# JRA55, on the PanArcticShenanigans grid. Built from the branch's examples/arctic_nested_ocean_sea_ice.jl with our
# domain, dates and forcing choices:
#   - the conformal conic grid (common.jl), Atlantic and Pacific edges open;
#   - GLORYS daily means and JRA55 for January 1993, as in our runs;
#   - JRA55 rivers (JRA55PrescribedLand), as in panarctic.jl;
#   - sea-level pressure set to 101325 Pa everywhere unless PRESSURE=true: GLORYS's η has no inverse-barometer
#     response, so this is the configuration to compare with GLORYS (see claude/pressure-forcing-glorys.md).
#
# Same environment as panarctic_nested_1deg.jl. Run:
#   SPACING_KM=100 DAYS=5 julia --project=. panarctic_nested_coupled.jl [output_dir]
# Environment: SPACING_KM (default 100), DAYS (default 5), PRESSURE (default false), DT_MINUTES (default 20 at 100 km;
# use 10 at 15 km).
#
# Compare nested_monitor.csv with stage 1 (ocean only) and, at 15 km, with our coupled-wind-nopressure-48h run.

include("common.jl")
using CopernicusMarine
using CUDA
using Dates

output_dir = isempty(ARGS) ? joinpath(@__DIR__, "output-coupled") : ARGS[1]
mkpath(output_dir)

spacing  = parse(Float64, get(ENV, "SPACING_KM", "100")) * kilometers
days     = parse(Float64, get(ENV, "DAYS", "5"))
pressure = get(ENV, "PRESSURE", "false") == "true"
Δt       = parse(Float64, get(ENV, "DT_MINUTES", spacing > 50kilometers ? "20" : "10")) * minutes

arch = GPU()
grid = panarctic_grid(arch; spacing)

start_date = DateTime(1993, 1, 1)
end_date   = start_date + Day(ceil(Int, days) + 1)
dates      = (start_date, end_date)

# ## Ocean and sea ice, nested in GLORYS

parent_dataset = GLORYSDaily()
ocean   = nested_ocean_simulation(grid, parent_dataset; dates, open_sides = (:south, :north))
sea_ice = nested_sea_ice_simulation(ocean, parent_dataset; dates)
child   = ocean.model.child

# ## JRA55 atmosphere, radiation and rivers
#
# With the pressure off, the whole record is held in memory so the constant pressure is not overwritten by a reload.

jra55 = (; dataset = MultiYearJRA55(), start_date, end_date)
atmosphere = JRA55PrescribedAtmosphere(arch; jra55...)
if !pressure
    atmosphere = JRA55PrescribedAtmosphere(arch; jra55..., time_indices_in_memory = length(atmosphere.velocities.u.times))
    fill!(parent(atmosphere.pressure.data), 101325)
    @warn "JRA55 sea-level pressure set to 101325 Pa everywhere (PRESSURE=false)"
end
radiation = JRA55PrescribedRadiation(arch; jra55...)
land      = JRA55PrescribedLand(child.grid; jra55...)

coupled_model = EarthSystemModel(; ocean, sea_ice, land, atmosphere, radiation)
simulation    = Simulation(coupled_model; Δt, stop_time = days * 24hours)

add_nested_monitor!(simulation, child, output_dir; exchanger = ocean.model.exchanger,
                    parent_η = parent_free_surface(grid, parent_dataset, dates))
add_nested_output!(simulation, child, output_dir)

function ice_progress(sim)
    h = sea_ice.model.ice_thickness
    ℵ = sea_ice.model.ice_concentration
    uᵢ, vᵢ = sea_ice.model.velocities
    @info @sprintf("         ice: max h %.2f m, mean ℵ %.3f, max |uᵢ|, |vᵢ| %.2f, %.2f m/s",
                   maximum(h), sum(interior(ℵ)) / length(interior(ℵ)), maximum(abs, uᵢ), maximum(abs, vᵢ))
    return nothing
end
add_callback!(simulation, ice_progress, TimeInterval(6hours))

run!(simulation)
