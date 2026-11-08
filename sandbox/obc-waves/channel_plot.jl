# Plot results/channel.jls (from channel_waves.jl) in the plotting environment:
#   julia --project=$OCEAN_HPC_ROOT/viz channel_plot.jl
using CairoMakie, Serialization
const g, H, f = 9.80665, 500.0, 1.4e-4
const c = sqrt(g * H)
theory_r(ω, α) = (c * ω / sqrt(ω^2 - f^2) - α) / (c * ω / sqrt(ω^2 - f^2) + α)

R = deserialize("results/channel.jls")
fig = Figure(size = (1300, 480), fontsize = 15)
ax = Axis(fig[1, 1], xlabel = "ω / f", ylabel = "reflection coefficient |r|",
          title = "(a) one Poincaré packet: Oceananigans against theory")
w = range(1.02, 8, length = 300)
lines!(ax, w, abs.(theory_r.(w .* f, c)), color = :steelblue, linewidth = 2, label = "Flather, theory")
hlines!(ax, [1.0], color = :firebrick, linewidth = 2, label = "Dirichlet, theory")
scatter!(ax, [p.ω / f for p in R[:packets]], [p.flather for p in R[:packets]], color = :steelblue, markersize = 14,
         label = "Flather, Oceananigans")
scatter!(ax, [p.ω / f for p in R[:packets]], [p.dirichlet for p in R[:packets]], color = :firebrick, markersize = 14,
         marker = :rect, label = "Dirichlet, Oceananigans")
ylims!(ax, 0, 1.05); axislegend(ax, position = :rc)

ax2 = Axis(fig[1, 2], xlabel = "distance from the open edge (km)", ylabel = "η (cm)",
           title = "(b) steady balanced jet, external η 2 cm high, day 5")
s = R[:flather_2cm]
lines!(ax2, s.y ./ 1e3, 100 .* s.truth, color = :black, linewidth = 3, label = "truth")
lines!(ax2, R[:dirichlet_2cm].y ./ 1e3, 100 .* R[:dirichlet_2cm].η, color = :firebrick, linestyle = :dash,
       linewidth = 2, label = "Dirichlet")
lines!(ax2, s.y ./ 1e3, 100 .* s.η, color = :steelblue, linewidth = 2, label = "Flather")
xlims!(ax2, 0, 2500); axislegend(ax2)
mkpath("figures"); save("figures/channel_validation.png", fig)
println("saved figures/channel_validation.png")
