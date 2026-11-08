# E2c: the external datum drifts slowly (3 cm over 10 days, run for 20 days), the realistic case that E2b (3 cm in
# one day) exaggerated. The question: does M6 (Flather, with C a slow average of the edge mismatch, so the net
# transport follows the external one at low frequency) stay exact on the datum while letting waves out as Flather does?
# M6 with averaging times of 12 h, 2 d and 5 d, against Flather, M5b and M7, on E2c and on the wave test E1.
#
#   julia -t 8 --project=$OCEAN_HPC_ROOT run_slow.jl                     # about 15 minutes
#   julia --project=$OCEAN_HPC_ROOT/viz compare.jl slow                  # figures/*_slow.png and the tables

include("obc_vortex.jl")
P0 = default_parameters()
isfile("results/reference_adjustment.jls") || run_reference(:adjustment; verbose = false)
for experiment in (:offset_slow, :adjustment)
    for m in (:flather, :balanced_net, :rotating_flather)
        experiment === :adjustment && m !== :flather && isfile("results/$(m)_adjustment.jls") && continue
        run_case(m, experiment; verbose = false)
        println("done $m $experiment")
    end
    for (suffix, T_C) in (("_T12h", 12hours), ("", 2days), ("_T5d", 5days))
        run_case(:balanced_filtered, experiment; P = merge(P0, (; T_C)), suffix, verbose = false)
        println("done balanced_filtered$suffix $experiment")
    end
end
println("Now plot with:  julia --project=\$OCEAN_HPC_ROOT/viz compare.jl slow")
