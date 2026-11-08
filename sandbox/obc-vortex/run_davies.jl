# Dirichlet with a sponge (M3, Davies) for different rim widths and relaxation times, against Flather, plain
# Dirichlet and M5b, on the wave test (E1, with its reference run) and the growing datum offset (E2b).
#
# On oceananigans-hpc (Oceananigans 0.117) the simulations and the plots use different environments:
#
#   julia -t 8 --project=$OCEAN_HPC_ROOT run_davies.jl                  # runs, results/*.jls
#   julia --project=$OCEAN_HPC_ROOT/viz compare.jl davies                # figures/*_davies.png and tables
#
# A barotropic wave near f has wavelength 2πc/f ≈ 3000 km here, so a 200 km rim cannot absorb it; the question
# is whether a much wider or gentler rim gets close to Flather on E1 while staying exact on E2b. The rim relaxes
# u and v only (the free surface takes no forcing), unlike a classic Davies zone, which also relaxes η.

include("obc_vortex.jl")
P0 = default_parameters()

rims = [("_R200",     (rim_width = 200kilometers,  rim_rate = 1 / 2hours)),
        ("_R500",     (rim_width = 500kilometers,  rim_rate = 1 / 2hours)),
        ("_R1000",    (rim_width = 1000kilometers, rim_rate = 1 / 2hours)),
        ("_R1000_6h", (rim_width = 1000kilometers, rim_rate = 1 / 6hours))]

for experiment in (:adjustment, :offset_ramp)
    experiment === :adjustment && run_reference(experiment; verbose = false)
    for m in (:flather, :dirichlet, :balanced_net)
        run_case(m, experiment; verbose = false)
        println("done $m $experiment")
    end
    for (suffix, rim) in rims
        run_case(:davies, experiment; P = merge(P0, rim), suffix, verbose = false)
        println("done davies$suffix $experiment")
    end
end

println("Now plot with:  julia --project=\$OCEAN_HPC_ROOT/viz compare.jl davies")
