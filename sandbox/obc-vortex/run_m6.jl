# M6 (filtered C) on E0, E2, E1 with T_C = 12 h, 2 days and 5 days, then the comparisons.
include("obc_vortex.jl")
P0 = default_parameters()
for experiment in (:steady, :offset, :adjustment)
    run_case(:balanced_filtered, experiment; verbose = false)
    experiment === :steady && continue
    run_case(:balanced_filtered, experiment; P = merge(P0, (T_C = 12hours,)), suffix = "_T12h", verbose = false)
    run_case(:balanced_filtered, experiment; P = merge(P0, (T_C = 5days,)),   suffix = "_T5d",  verbose = false)
end
include("analysis.jl")
for experiment in (:steady, :offset, :adjustment)
    compare(experiment)
end
for experiment in (:offset, :adjustment)
    compare(experiment; methods = ("flather", "balanced_net", "balanced_filtered_T12h", "balanced_filtered", "balanced_filtered_T5d"),
            tag_suffix = "_M6_sensitivity")
end
