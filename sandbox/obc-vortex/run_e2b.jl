# E2b: a 3 cm datum offset that grows over one day (the case that defeats the time-filtered C).
include("obc_vortex.jl")
P0 = default_parameters()
for m in (:flather, :dirichlet, :balanced_net, :balanced_pv)
    run_case(m, :offset_ramp; verbose = false)
end
run_case(:balanced_filtered, :offset_ramp; verbose = false)
run_case(:balanced_filtered, :offset_ramp; P = merge(P0, (T_C = 12hours,)), suffix = "_T12h", verbose = false)
run_case(:balanced_filtered, :offset_ramp; P = merge(P0, (T_C = 5days,)),   suffix = "_T5d",  verbose = false)
include("analysis.jl")
compare(:offset_ramp; methods = ("flather", "dirichlet", "balanced_net", "balanced_pv",
                                 "balanced_filtered_T12h", "balanced_filtered", "balanced_filtered_T5d"))
