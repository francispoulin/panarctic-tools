# Transport-trusting methods (and Flather as the baseline) on the shelf (E4) and the slow net inflow (E3).
#   julia --project=... run_trust.jl shelf        # E1, E2, E2b with the shelf, plus the shelf reference
#   julia --project=... run_trust.jl net_inflow   # E3, 30 days, 6-hourly output
include("obc_vortex.jl")
methods = (:flather, :dirichlet, :davies, :flather_target, :balanced_net, :balanced_pv)
part = isempty(ARGS) ? "shelf" : ARGS[1]
if part == "shelf"
    run_reference(:adjustment; shelf = true, verbose = false)
    for e in (:adjustment, :offset, :offset_ramp), m in methods
        run_case(m, e; shelf = true, verbose = false); println("done $m $e shelf")
    end
    include("analysis.jl")
    for e in (:adjustment, :offset, :offset_ramp)
        compare(e; shelf = true, methods = String.(methods))
    end
else
    for m in methods
        run_case(m, :net_inflow; output_interval = 6hours, verbose = false); println("done $m net_inflow")
    end
    include("analysis.jl")
    compare(:net_inflow; methods = String.(methods))
end
