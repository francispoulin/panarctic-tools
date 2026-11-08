# M7, Flather with the rotational correction (see ../obc-waves), on the wave test (E1) and the growing datum (E2b),
# next to Flather and M5b. Run after run_davies.jl (it reuses the E1 reference), then summary.jl:
#
#   julia -t 8 --project=$OCEAN_HPC_ROOT run_rotating.jl
#   julia --project=$OCEAN_HPC_ROOT/viz compare.jl rotating
#   julia --project=$OCEAN_HPC_ROOT/viz summary.jl

include("obc_vortex.jl")
for experiment in (:adjustment, :offset_ramp)
    isfile("results/reference_adjustment.jls") || experiment !== :adjustment || run_reference(:adjustment; verbose = false)
    run_case(:rotating_flather, experiment; verbose = false)
    println("done rotating_flather $experiment")
end
println("Now plot with:  julia --project=\$OCEAN_HPC_ROOT/viz compare.jl rotating")
