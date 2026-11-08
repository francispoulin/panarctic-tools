# Figures and tables from results/*.jls, in a plotting environment (CairoMakie; no Oceananigans needed):
#
#   julia --project=$OCEAN_HPC_ROOT/viz compare.jl davies              # the rim-width test (run_davies.jl)
#   julia --project=$OCEAN_HPC_ROOT/viz compare.jl adjustment offset   # any experiments, all default methods

include("analysis.jl")

const davies_methods = ("flather", "dirichlet", "balanced_net",
                        "davies_R200", "davies_R500", "davies_R1000", "davies_R1000_6h")

args = isempty(ARGS) ? ["davies"] : ARGS
for a in args
    if a == "davies"
        for experiment in (:adjustment, :offset_ramp)
            compare(experiment; methods = davies_methods, tag_suffix = "_davies")
        end
    else
        compare(Symbol(a))
    end
end
