# Is the extra inflow at the Atlantic edge the Flather term? And where is the model's free surface off GLORYS?
#
# Include in a live session after (part of) a run, e.g. in the run1d REPL after run!(simulation):
#
#   include("../panarctic-tools/panarctic-edge-eta.jl")
#
# At the south (Atlantic) edge the gravity-wave radiation (Flather) condition in Oceananigans 0.113.6 sets the
# barotropic transport per unit width on the boundary faces to
#
#     Vᵇ = Vᵉˣᵗ − √(g H) (ηᵇ − ηᵉˣᵗ),      ηᵇ = (η[i, 0] + η[i, 1]) / 2,
#
# so the model's extra inflow relative to GLORYS should be −√(g H) (ηᵇ − ηᵉˣᵗ), column by column. This script
# compares the model's boundary transport with GLORYS and with that prediction, and maps η_model − η_GLORYS.
# GLORYS fields are interpolated linearly in time between daily records, as the boundary condition does.
#
# Writes edge_flather_h<hour>.png and eta_diff_h<hour>.png into output_dir.

using Oceananigans.Architectures: on_architecture
using Oceananigans.Operators: Δxᶜᶠᶜ, Δzᶜᶠᶜ

function edge_glorys_at(fts, t)
    ts = fts.times
    n  = clamp(searchsortedlast(ts, t), 1, length(ts) - 1)
    α  = clamp((t - ts[n]) / (ts[n+1] - ts[n]), 0, 1)
    a  = replace(Array(interior(fts[n])),   NaN => 0.0)
    b  = replace(Array(interior(fts[n+1])), NaN => 0.0)
    return (1 - α) .* a .+ α .* b
end

function edge_flather_report()
    t     = ocean.model.clock.time
    hour  = round(Int, t / hours)
    g     = ocean.model.free_surface.gravitational_acceleration
    ηc    = on_architecture(CPU(), ocean.model.free_surface.displacement).data   # OffsetArray, halos included
    kη    = size(ηc, 3) == 1 ? first(axes(ηc, 3)) : Nz + 1                     # η lives at the top face
    Vm    = Array(interior(ocean.model.free_surface.barotropic_velocities.V))[:, 1, 1]   # m²/s at the south faces
    vg    = edge_glorys_at(vᶠ, t)
    ηext  = edge_glorys_at(ηᵉˣᵗ, t)[:, 1, 1]

    Hs  = [sum(wetcell(i, 1, k, cpu_grid, Center(), Face(), Center()) ? Δzᶜᶠᶜ(i, 1, k, cpu_grid) : 0.0 for k in 1:Nz) for i in 1:Nx]
    Vg  = [sum(wetcell(i, 1, k, cpu_grid, Center(), Face(), Center()) ? vg[i, 1, k] * Δzᶜᶠᶜ(i, 1, k, cpu_grid) : 0.0 for k in 1:Nz) for i in 1:Nx]
    Δx  = [Δxᶜᶠᶜ(i, 1, Nz, cpu_grid) for i in 1:Nx]
    wet = Hs .> 0
    ηb  = [(ηc[i, 0, kη] + ηc[i, 1, kη]) / 2 for i in 1:Nx]
    c   = sqrt.(g .* max.(Hs .+ ηb, 0))
    Δη  = ifelse.(wet, ηb .- ηext, 0.0)
    flather = -c .* Δη                                  # predicted model − GLORYS, m²/s

    Sv(a) = sum(a[wet] .* Δx[wet]) / 1e6
    width = sum(Δx[wet])
    @printf("\nAtlantic edge at hour %d (%d wet columns, %.0f km)\n", hour, count(wet), width / 1e3)
    @printf("  net in: model %+7.2f Sv, GLORYS %+7.2f Sv, difference %+7.2f Sv\n", Sv(Vm), Sv(Vg), Sv(Vm .- Vg))
    @printf("  Flather term −√(gH)(ηᵇ − ηᵉˣᵗ) summed: %+7.2f Sv   (should match the difference)\n", Sv(flather))
    @printf("  edge-mean ηᵇ − ηᵉˣᵗ: %+.2f cm (width-weighted); range %+.1f to %+.1f cm\n",
            100 * sum(Δη[wet] .* Δx[wet]) / width, 100 * minimum(Δη[wet]), 100 * maximum(Δη[wet]))
    for (name, r) in (("i < 250", 1:249), ("250 ≤ i < 335", 250:334), ("335 ≤ i ≤ 352", 335:352), ("i > 352", 353:Nx))
        m = falses(Nx); m[r] .= true; m .&= wet
        any(m) || continue
        @printf("  %-14s  model − GLORYS %+6.2f Sv, Flather term %+6.2f Sv, mean Δη %+5.2f cm\n", name,
                sum((Vm .- Vg)[m] .* Δx[m]) / 1e6, sum(flather[m] .* Δx[m]) / 1e6, 100 * sum(Δη[m] .* Δx[m]) / sum(Δx[m]))
    end

    is = findall(wet)
    fig = Figure(size = (1100, 750))
    ax = Axis(fig[1, 1], ylabel = "model − GLORYS into domain\nper column [Sv]", title = "Atlantic edge, hour $hour")
    lines!(ax, is, ((Vm .- Vg) .* Δx ./ 1e6)[is], label = "model − GLORYS")
    lines!(ax, is, (flather .* Δx ./ 1e6)[is], label = "−√(gH)(ηᵇ − ηᵉˣᵗ)", linestyle = :dash)
    axislegend(ax, position = :lt)
    ax2 = Axis(fig[2, 1], xlabel = "i", ylabel = "η [m]")
    lines!(ax2, is, ηb[is], label = "model ηᵇ")
    lines!(ax2, is, ηext[is], label = "GLORYS ηᵉˣᵗ")
    axislegend(ax2, position = :lt)
    save(joinpath(output_dir, "edge_flather_h$(hour).png"), fig)

    # Whole-domain map of η_model − η_GLORYS
    ηm  = Array(interior(ocean.model.free_surface.displacement))[:, :, 1]
    ηg  = edge_glorys_at(ηᵉˣᵗ, t)[:, :, 1]
    top = [wetcell(i, j, Nz, cpu_grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
    d   = ifelse.(top, ηm .- ηg, NaN)
    fig = Figure(size = (900, 800))
    ax = Axis(fig[1, 1], title = "η model − GLORYS [m], hour $hour", xlabel = "i", ylabel = "j", aspect = DataAspect())
    hm = heatmap!(ax, d, colormap = :balance, colorrange = (-0.15, 0.15), nan_color = :gray80)
    Colorbar(fig[1, 2], hm)
    save(joinpath(output_dir, "eta_diff_h$(hour).png"), fig)
    println("Saved edge_flather_h$(hour).png and eta_diff_h$(hour).png in $output_dir")
    return nothing
end

edge_flather_report()
