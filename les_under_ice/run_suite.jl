# Runs several cases one after another: the LES (skipped if its output exists), the CATKE and k-ε columns, the NEMO TKE
# columns (nemo_tke.jl, seconds on the CPU), and the comparisons
# (compare.jl against our LES; compare_lesbrary.jl against the calibration LES for the _stokes cases). In a REPL on a
# GPU, from the PanArcticShenanigans clone's environment:
#
#     empty!(ARGS); push!(ARGS, "med_wind_med_cooling", "med_wind_med_cooling_stokes")      # or no arguments: all new
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/run_suite.jl")
#
# No arguments: the three wind-and-cooling cases of the 24-hour suite, without and with waves (6 LES, ~50 min each).
# A case whose LES output exists is not rerun (delete output/<case>_les.jld2 to rerun it).

include(joinpath(@__DIR__, "cases.jl"))

suite_cases = isempty(ARGS) ? ["weak_wind_strong_cooling", "med_wind_med_cooling", "strong_wind_weak_cooling",
                               "weak_wind_strong_cooling_stokes", "med_wind_med_cooling_stokes",
                               "strong_wind_weak_cooling_stokes"] : copy(ARGS)

function run_script(script, args...)
    empty!(ARGS)
    append!(ARGS, args)
    Base.include(Main, joinpath(@__DIR__, script))
    return nothing
end

for name in suite_cases
    haskey(cases, name) || (@warn "Unknown case $name: skipped"; continue)
    if isfile(joinpath(output_folder(), name * "_les.jld2"))
        @info "LES $name exists: not rerun"
    else
        run_script("les.jl", name)
    end
    run_script("columns.jl", name)
    run_script("nemo_tke.jl", name)
    run_script("compare.jl", name)
    endswith(name, "_stokes") && run_script("compare_lesbrary.jl", replace(name, "_stokes" => ""))
end
empty!(ARGS)
@info "Suite done: $(join(suite_cases, ", ")). Tables in $(output_folder())"
