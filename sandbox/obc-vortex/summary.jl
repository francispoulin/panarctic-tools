# One table and one figure that tell the rim-width story (run after run_davies.jl):
#
#   julia --project=$OCEAN_HPC_ROOT/viz summary.jl            # figures/obc_tradeoff.png, table on screen
#
# Every condition is scored on the two failures that matter, each against the known truth at day 2:
#   datum error  (E2b): the external η is 2 cm too high but its transport is right; the truth is the steady vortex.
#                       Score: rms η error, and the peak net-transport error during the run.
#   wave error   (E1):  an unbalanced bump adjusts and its waves must leave; the truth is the extended-basin run.
#                       The η error is split into its basin-mean part (volume that should have left or stayed out)
#                       and the rest (the shape: waves trapped in the basin), and the unbalanced energy E′ is given
#                       as a multiple of the truth's own E′ (1 = waves leave as in the truth).

include("analysis.jl")

const summary_methods = ["flather", "flather_sponge_R200", "flather_sponge_R500",
                         "dirichlet", "davies_R200", "davies_R500", "davies_R1000", "davies_R1000_6h",
                         "balanced_net", "balanced_net_sponge_R200", "balanced_net_sponge_R500", "rotating_flather"]

function scores(name, experiment)
    r = load(name, experiment)
    P = r.P
    refs = reference_states(experiment, P, r.times)
    n = length(r.snapshots)
    s, ref = r.snapshots[n], refs[min(n, length(refs))]
    Δη = Float64.(s.η) .- ref.η
    m = sum(Δη) / length(Δη)
    inv = basin_inverter(P)
    E′    = split_energy(inv, Float64.(s.u), Float64.(s.v), Float64.(s.η); v_ext = Float64.(s.v[:, 1])).E′
    E′ref = split_energy(inv, ref.u, ref.v, ref.η; v_ext = Float64.(ref.v[:, 1])).E′
    return (total = 100rms(Δη), mean = 100m, shape = 100rms(Δη .- m), E′ratio = E′ / E′ref,
            peak = maximum(abs, r.net .- r.net_ext) / 1e6, mean_true = 100sum(ref.η) / length(ref.η))
end

rows = [(m, scores(m, :adjustment), scores(m, :offset_ramp)) for m in summary_methods]

println("\nTruth: basin-mean η at day 2 = $(round(rows[1][2].mean_true, digits = 3)) cm (E1), ",
        "$(round(rows[1][3].mean_true, digits = 3)) cm (E2b)\n")
@printf("| %-26s | %8s | %8s | %8s | %8s | %9s | %9s |\n", "method", "E2b rms", "E2b peak",
        "E1 rms", "E1 mean", "E1 shape", "E1 E′/E′true")
@printf("| %-26s | %8s | %8s | %8s | %8s | %9s | %9s |\n", "", "cm", "Sv", "cm", "cm", "cm", "")
println("|", "-"^28, ("|" * "-"^10)^4, ("|" * "-"^11)^2, "|")
for (m, a, o) in rows
    @printf("| %-26s | %8.2f | %8.2f | %8.2f | %+8.2f | %9.2f | %9.1f |\n",
            label_of(m), o.total, o.peak, a.total, a.mean, a.shape, a.E′ratio)
end

# Figure: (a) the trade-off, one point per condition; (b) the E1 error split into volume and shape.
fig = Figure(size = (1400, 560), fontsize = 15)
ax = Axis(fig[1, 1], xlabel = "wave error: E1 rms η error, shape part (cm)",
          ylabel = "datum error: E2b rms η error (cm)",
          title = "(a) each condition against the two failures (best = bottom left)")
for (m, a, o) in rows
    scatter!(ax, [a.shape], [o.total]; color = colour_of(m), markersize = 16)
    text!(ax, a.shape, o.total; text = "  " * label_of(m), fontsize = 11, align = (:left, :center))
end
xlims!(ax, 0, 0.6); ylims!(ax, -0.15, 2.6)

ax2 = Axis(fig[1, 2], ylabel = "E1 rms η error (cm)", title = "(b) E1 error: volume part and shape part",
           xticks = (1:length(rows), [label_of(m) for (m, _, _) in rows]), xticklabelrotation = π / 3)
xs = 1:length(rows)
barplot!(ax2, xs, [abs(a.mean) for (_, a, _) in rows]; color = :gray70, label = "basin-mean part |⟨Δη⟩|")
barplot!(ax2, xs, [a.shape for (_, a, _) in rows]; offset = [abs(a.mean) for (_, a, _) in rows],
         color = [colour_of(m) for (m, _, _) in rows], label = "shape part (trapped waves)")
axislegend(ax2, position = :lt)

mkpath("figures")
save("figures/obc_tradeoff.png", fig)
println("\n  saved figures/obc_tradeoff.png")
println("  (the stacked bar adds the two parts; the rms of the total is √(mean² + shape²), slightly less)")
