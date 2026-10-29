# Local diagnostics, not for committing. Paste into panarctic.jl just before `run!(simulation)`.
# Every 6 hours: where are the fastest currents? Daily: full 3-D ocean output for the strait and window tools.

function report_maxima(sim)
    u, v, w = sim.model.ocean.model.velocities
    ℓ = (Center(), Center(), Center())
    for (name, φ) in (("max|u|", u), ("max|v|", v))
        a = abs.(Array(interior(φ)))
        replace!(a, NaN => 0)
        m, I = findmax(a)
        i, j, k = Tuple(I)
        iᶜ, jᶜ = min(i, Nx), min(j, Ny)
        @info @sprintf("    %s = %.2f m/s at (%d, %d, %d) ≈ %.1f°E %.1f°N, z = %.0f m, %d cells from an open edge",
                       name, m, i, j, k,
                       Oceananigans.Grids.λnode(iᶜ, jᶜ, k, cpu_grid, ℓ...), Oceananigans.Grids.φnode(iᶜ, jᶜ, k, cpu_grid, ℓ...),
                       Oceananigans.Grids.znode(iᶜ, jᶜ, k, cpu_grid, ℓ...), min(j - 1, Ny - j))
    end
    return nothing
end

add_callback!(simulation, report_maxima, TimeInterval(6hours))

ocean.output_writers[:full] = JLD2Writer(ocean.model, merge(ocean.model.tracers, (; u, v, w));
                                         filename = joinpath(output_dir, "arctic_ocean_full.jld2"),
                                         schedule = TimeInterval(1days),
                                         overwrite_files = true)
