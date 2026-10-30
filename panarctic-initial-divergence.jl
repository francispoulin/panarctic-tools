# Is the initial velocity field consistent with the initial free surface on the model's grid?
#
# The model starts from GLORYS u, v rotated and interpolated onto its own grid and bathymetry. On GLORYS's grid
# those velocities satisfy ∂η/∂t = −∇·(∫u dz); on ours they need not. This script computes, at t = 0,
#
#   (a) the free-surface tendency implied by the mapped velocities:  ∂η/∂t = −∇·U  (U = depth-integrated transport
#       through the model's own wet faces), and
#   (b) GLORYS's own tendency between its first two daily records, on the model grid,
#
# and maps both and their difference (the spurious part). No model state is needed: include it in any session that
# has run panarctic_setup.jl (uses uᶠ, vᶠ, ηᵉˣᵗ, cpu_grid, wetcell, Nx, Ny, Nz, output_dir):
#
#   include("../panarctic-tools/panarctic-initial-divergence.jl")
#
# Writes initial_divergence.png into output_dir.

using Oceananigans.Operators: Δxᶜᶠᶜ, Δzᶜᶠᶜ, Δyᶠᶜᶜ, Δzᶠᶜᶜ, Azᶜᶜᶜ

function initial_divergence_report(; n = 1)
    u0 = replace(Array(interior(uᶠ[n])), NaN => 0.0)     # (Nx+1, Ny, Nz) in the grid frame
    v0 = replace(Array(interior(vᶠ[n])), NaN => 0.0)     # (Nx, Ny+1, Nz)
    Nxu, Nyv = size(u0, 1), size(v0, 2)

    Fu = zeros(Nxu, Ny)                                  # m³/s through each u-face column
    for j in 1:Ny, i in 1:Nxu, k in 1:Nz
        wetcell(i, j, k, cpu_grid, Face(), Center(), Center()) || continue
        Fu[i, j] += u0[i, j, k] * Δyᶠᶜᶜ(i, j, k, cpu_grid) * Δzᶠᶜᶜ(i, j, k, cpu_grid)
    end
    Fv = zeros(Nx, Nyv)                                  # m³/s through each v-face column
    for j in 1:Nyv, i in 1:Nx, k in 1:Nz
        wetcell(i, j, k, cpu_grid, Center(), Face(), Center()) || continue
        Fv[i, j] += v0[i, j, k] * Δxᶜᶠᶜ(i, j, k, cpu_grid) * Δzᶜᶠᶜ(i, j, k, cpu_grid)
    end

    top = [wetcell(i, j, Nz, cpu_grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
    Az  = [Azᶜᶜᶜ(i, j, Nz, cpu_grid) for i in 1:Nx, j in 1:Ny]
    cmday = 100 * 86400

    implied = [top[i, j] ? -(Fu[min(i + 1, Nxu), j] - Fu[i, j] + Fv[i, min(j + 1, Nyv)] - Fv[i, j]) / Az[i, j] : NaN
               for i in 1:Nx, j in 1:Ny] .* cmday

    η1 = replace(Array(interior(ηᵉˣᵗ[n]))[:, :, 1],     NaN => 0.0)
    η2 = replace(Array(interior(ηᵉˣᵗ[n + 1]))[:, :, 1], NaN => 0.0)
    actual = ifelse.(top, (η2 .- η1) ./ (ηᵉˣᵗ.times[n + 1] - ηᵉˣᵗ.times[n]) .* cmday, NaN)
    spurious = implied .- actual

    w = top .* Az
    mean_(a) = sum(ifelse.(top, a, 0.0) .* w) / sum(w)
    rms_(a)  = sqrt(mean_(ifelse.(top, a, 0.0) .^ 2))
    @printf("\nFree-surface tendency at t = 0 [cm/day], area-weighted over %d wet columns:\n", count(top))
    @printf("  implied by mapped velocities:  mean %+7.2f   rms %7.2f\n", mean_(implied), rms_(implied))
    @printf("  GLORYS records %d → %d:          mean %+7.2f   rms %7.2f\n", n, n + 1, mean_(actual), rms_(actual))
    @printf("  difference (spurious):         mean %+7.2f   rms %7.2f\n", mean_(spurious), rms_(spurious))
    @printf("  net transport in through the open edges from the mapped velocities: %+6.2f Sv\n",
            (sum(Fv[:, 1]) - sum(Fv[:, Nyv])) / 1e6)

    fig = Figure(size = (1800, 650))
    for (c, (title, a)) in enumerate((("∂η/∂t implied by mapped u, v at t = 0", implied),
                                      ("∂η/∂t of GLORYS, records $n → $(n + 1)", actual),
                                      ("difference (spurious)", spurious)))
        ax = Axis(fig[1, c], title = title, xlabel = "i", ylabel = "j", aspect = DataAspect())
        hm = heatmap!(ax, a, colormap = :balance, colorrange = (-50, 50), nan_color = :gray80)
        c == 3 && Colorbar(fig[1, 4], hm, label = "cm/day")
    end
    path = joinpath(output_dir, "initial_divergence.png")
    save(path, fig)
    println("Saved $path")
    return (; implied, actual, spurious, Fu, Fv)
end

initial_divergence = initial_divergence_report()
