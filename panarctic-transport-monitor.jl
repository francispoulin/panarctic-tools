# Hourly net transport through chosen rows of v-faces, model vs GLORYS, and the domain-mean free surface.
#
# Include AFTER panarctic_setup.jl (it uses grid, cpu_grid, wetcell, vᶠ, ηᵉˣᵗ, ocean, simulation, output_dir)
# and BEFORE run!(simulation):
#
#   include("panarctic_setup.jl")
#   include("../panarctic-tools/panarctic-transport-monitor.jl")
#   run!(simulation)
#
# GLORYS transports use GLORYS velocities (already rotated into the grid frame, on the model's faces) times the
# model's own face areas, interpolated linearly in time between the daily records. Because the model starts from
# vᶠ[1] and ηᵉˣᵗ[1], model and "GLORYS" agree at t = 0 by construction, so their difference afterwards is drift.
#
# Prints one line per hour and writes everything to <output_dir>/transport_monitor.csv.

using Oceananigans.Operators: Δxᶜᶠᶜ, Δzᶜᶠᶜ, Azᶜᶜᶜ
using Oceananigans.Architectures: on_architecture

# (j, sign): sign makes positive = into the domain at the open edges, and towards +j for interior rows.
monitor_rows   = [(1, +1), (70, +1), (130, +1), (175, +1), (186, +1), (Ny + 1, -1)]
monitor_window = 335:352      # Norwegian jet band, reported at j = 175
monitor_window_row = 175

monitor_row_area(j) = [wetcell(i, j, k, cpu_grid, Center(), Face(), Center()) ?
                       Δxᶜᶠᶜ(i, j, k, cpu_grid) * Δzᶜᶠᶜ(i, j, k, cpu_grid) : 0.0
                       for i in 1:Nx, _ in 1:1, k in 1:Nz]

monitor_A_cpu = [monitor_row_area(j) for (j, _) in monitor_rows]
monitor_W_cpu = let A = monitor_row_area(monitor_window_row)
    A[setdiff(1:Nx, monitor_window), :, :] .= 0
    A
end
monitor_Aη_cpu = [wetcell(i, j, Nz, cpu_grid, Center(), Center(), Center()) ? Azᶜᶜᶜ(i, j, Nz, cpu_grid) : 0.0
                  for i in 1:Nx, j in 1:Ny, _ in 1:1]

monitor_A  = [on_architecture(arch, A) for A in monitor_A_cpu]
monitor_W  = on_architecture(arch, monitor_W_cpu)
monitor_Aη = on_architecture(arch, monitor_Aη_cpu)
monitor_area_total = sum(monitor_Aη_cpu)

# GLORYS, once per daily record (CPU)
monitor_times = collect(vᶠ.times)
monitor_glorys_rows = zeros(length(monitor_rows), length(monitor_times))
monitor_glorys_win  = zeros(length(monitor_times))
monitor_glorys_η    = zeros(length(monitor_times))
for n in eachindex(monitor_times)
    local vg = replace(Array(interior(vᶠ[n])), NaN => 0.0)
    for (r, (j, s)) in enumerate(monitor_rows)
        monitor_glorys_rows[r, n] = s * sum(vg[:, j:j, :] .* monitor_A_cpu[r])
    end
    monitor_glorys_win[n] = sum(vg[:, monitor_window_row:monitor_window_row, :] .* monitor_W_cpu)
    local ηg = replace(Array(interior(ηᵉˣᵗ[n])), NaN => 0.0)
    monitor_glorys_η[n] = sum(ηg .* monitor_Aη_cpu) / monitor_area_total
end

function monitor_interp(series, t)
    ts = monitor_times
    t ≤ ts[1]   && return series[1]
    t ≥ ts[end] && return series[end]
    n = searchsortedlast(ts, t)
    α = (t - ts[n]) / (ts[n+1] - ts[n])
    return (1 - α) * series[n] + α * series[n+1]
end

monitor_file = joinpath(output_dir, "transport_monitor.csv")
open(monitor_file, "w") do io
    cols = ["hour"]
    for (j, _) in monitor_rows
        push!(cols, "j$(j)_model_Sv", "j$(j)_glorys_Sv")
    end
    push!(cols, "jet_model_Sv", "jet_glorys_Sv", "eta_model_m", "eta_glorys_m")
    println(io, join(cols, ","))
end

@printf("\n%6s  %15s  %15s  %15s  %15s  %17s\n", "hour", "Atlantic in", "Bering in", "j=175 net", "jet j=175",
        "mean η [m]")
@printf("%6s  %15s  %15s  %15s  %15s  %17s\n", "", "model / GLORYS", "model / GLORYS", "model / GLORYS",
        "model / GLORYS", "model / GLORYS")

function transport_monitor(sim)
    t   = ocean.model.clock.time
    vm  = interior(ocean.model.velocities.v)
    ηm  = interior(ocean.model.free_surface.displacement)

    model_rows  = [s * sum(vm[:, j:j, :] .* monitor_A[r]) for (r, (j, s)) in enumerate(monitor_rows)] ./ 1e6
    glorys_rows = [monitor_interp(view(monitor_glorys_rows, r, :), t) for r in eachindex(monitor_rows)] ./ 1e6
    model_win   = sum(vm[:, monitor_window_row:monitor_window_row, :] .* monitor_W) / 1e6
    glorys_win  = monitor_interp(monitor_glorys_win, t) / 1e6
    model_η     = sum(ηm .* monitor_Aη) / monitor_area_total
    glorys_η    = monitor_interp(monitor_glorys_η, t)

    r175 = findfirst(r -> r[1] == 175, monitor_rows)
    @printf("%6.1f  %+6.2f / %+6.2f  %+6.2f / %+6.2f  %+6.2f / %+6.2f  %+6.2f / %+6.2f  %+7.3f / %+7.3f\n",
            t / hours, model_rows[1], glorys_rows[1], model_rows[end], glorys_rows[end],
            model_rows[r175], glorys_rows[r175], model_win, glorys_win, model_η, glorys_η)

    open(monitor_file, "a") do io
        vals = [t / hours]
        for r in eachindex(monitor_rows)
            push!(vals, model_rows[r], glorys_rows[r])
        end
        push!(vals, model_win, glorys_win, model_η, glorys_η)
        println(io, join(vals, ","))
    end
    return nothing
end

transport_monitor(simulation)                                   # t = 0: model and GLORYS should agree
add_callback!(simulation, transport_monitor, TimeInterval(1hours))
@info "Transport monitor on: $(monitor_file)"
