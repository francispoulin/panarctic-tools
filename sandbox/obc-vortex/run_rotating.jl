# M7 (rotating Flather, see ../obc-waves) and M5b + M7 (the net constraint with the rotational term), next to
# Flather and M5b, on the wave test (E1) and the datum tests (E2b: 3 cm over 1 day; E2c: over 10 days).
#
#   julia -t 8 --project=$OCEAN_HPC_ROOT run_rotating.jl
#   julia --project=$OCEAN_HPC_ROOT/viz compare.jl rotating

include("obc_vortex.jl")
isfile("results/reference_adjustment.jls") || run_reference(:adjustment; verbose = false)
for experiment in (:adjustment, :offset_ramp, :offset_slow), m in (:flather, :balanced_net, :rotating_flather, :balanced_rotating)
    isfile("results/$(m)_$(experiment).jls") && m !== :balanced_rotating && continue
    run_case(m, experiment; verbose = false)
    println("done $m $experiment")
end
println("Now plot with:  julia --project=\$OCEAN_HPC_ROOT/viz compare.jl rotating")
