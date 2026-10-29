# Local debugging, not for committing. Paste into panarctic.jl just before `run!(simulation)`.
# Takes the first 6 time steps one at a time and reports, after each, which fields hold NaN and where the
# first one is; then stops the script, so the full run does not start.

function report_nans(sim)
    ocean = sim.model.ocean.model
    ice   = sim.model.sea_ice.model
    fields = (T  = ocean.tracers.T,           S  = ocean.tracers.S,
              u  = ocean.velocities.u,        v  = ocean.velocities.v,
              η  = ocean.free_surface.η,
              h  = ice.ice_thickness,         ℵ  = ice.ice_concentration,
              uᵢ = ice.velocities.u,          vᵢ = ice.velocities.v)
    ℓ = (Center(), Center(), Center())
    @info "After iteration $(iteration(sim)):"
    for (name, f) in pairs(fields)
        a = Array(interior(f))
        n = count(isnan, a)
        if n == 0
            @info "    $name: no NaN"
        else
            I = Tuple(findfirst(isnan, a))
            i, j = min(I[1], Nx), min(I[2], Ny)
            @info @sprintf("    %s: %d NaN, first at %s ≈ %.1f°E %.1f°N", name, n, string(I),
                           Oceananigans.Grids.λnode(i, j, 1, cpu_grid, ℓ...), Oceananigans.Grids.φnode(i, j, 1, cpu_grid, ℓ...))
        end
    end
    return nothing
end

report_nans(simulation)
for n in 1:6
    simulation.stop_iteration = n
    run!(simulation)
    report_nans(simulation)
end
exit()
