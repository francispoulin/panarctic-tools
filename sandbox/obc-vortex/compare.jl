# Figures and tables from results/*.jls, in a plotting environment (CairoMakie; no Oceananigans needed):
#
#   julia --project=$OCEAN_HPC_ROOT/viz compare.jl davies              # the rim-width test (run_davies.jl)
#   julia --project=$OCEAN_HPC_ROOT/viz compare.jl adjustment offset   # any experiments, all default methods

include("analysis.jl")

const davies_methods = ("flather", "dirichlet", "balanced_net",
                        "davies_R200", "davies_R500", "davies_R1000", "davies_R1000_6h",
                        "flather_sponge_R200", "flather_sponge_R500",
                        "balanced_net_sponge_R200", "balanced_net_sponge_R500")

const slow_methods = ("flather", "balanced_net", "rotating_flather",
                      "balanced_filtered_T12h", "balanced_filtered", "balanced_filtered_T5d")

args = isempty(ARGS) ? ["davies"] : ARGS
for a in args
    if a == "davies"
        for experiment in (:adjustment, :offset_ramp)
            compare(experiment; methods = davies_methods, tag_suffix = "_davies")
        end
    elseif a == "slow"
        compare(:offset_slow; methods = slow_methods, tag_suffix = "_slow")
        compare(:adjustment; methods = slow_methods, tag_suffix = "_slow")
    elseif a == "rotating"
        for experiment in (:adjustment, :offset_ramp, :offset_slow)
            compare(experiment; methods = ("flather", "balanced_net", "rotating_flather", "balanced_rotating"),
                    tag_suffix = "_rotating")
        end
    else
        compare(Symbol(a))
    end
end
