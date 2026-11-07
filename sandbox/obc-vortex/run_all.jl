# Run every method on the chosen experiments, plus the reference run where one is needed, then the comparison.
#
#   julia --project=<env with Oceananigans and CairoMakie> run_all.jl                      # E0, E2, E1
#   julia --project=... run_all.jl adjustment offset net_inflow                         # chosen experiments
#   julia --project=... run_all.jl --shelf adjustment offset                            # E4: with the shelf
#
# Output: results/*.jls, figures/*.png

include("obc_vortex.jl")

shelf = "--shelf" in ARGS
names = filter(a -> !startswith(a, "--"), ARGS)
experiments = isempty(names) ? (:steady, :offset, :adjustment) : Symbol.(names)

for experiment in experiments
    println("\n=== $experiment$(shelf ? " (shelf)" : "") ===")
    experiment === :adjustment && run_reference(experiment; shelf)
    for method in all_methods
        run_case(method, experiment; shelf, verbose = false)
        println("  done $method")
    end
end

include("analysis.jl")
for experiment in experiments
    compare(experiment; shelf)
end
