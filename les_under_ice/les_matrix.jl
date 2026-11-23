# Runs the calibration LES matrix (cases.jl, `les_matrix`: 21 cases, 1 m, 4 days; rotating stress 1 day) on one GPU.
# Start it in up to three REPLs, one per GPU (pelee and rondeau share the folder): each takes the next case nobody has
# claimed, runs it, and moves on, so the GPUs share the work however long each case takes.
#
#     empty!(ARGS); push!(ARGS, "2")          # the GPU, numbered as in nvidia-smi
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/les_matrix.jl")
#
# Bookkeeping in output/les_matrix/: <case>.claimed (host, GPU, start time), <case>.done (end time) or <case>.failed
# (the error). A finished case is skipped; a failed one is not retried (delete its .failed and .claimed to retry).
# A run stopped by Ctrl-C or a reboot leaves a .claimed without .done: delete that .claimed to run it again.
# Each LES runs in a fresh module (les.jl and cases.jl redefine globals and constants), and the GPU is freed after.
# Status of the whole matrix from any shell:
#     ls /mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/output/les_matrix/

using CUDA, Dates, Printf

matrix_gpu = isempty(ARGS) ? 0 : parse(Int, ARGS[1])
empty!(ARGS)
CUDA.device!(matrix_gpu)

matrix_folder = @__DIR__
include(joinpath(matrix_folder, "cases.jl"))                 # for the list les_matrix and output_folder()
matrix_state = mkpath(joinpath(output_folder(), "les_matrix"))
matrix_host = gethostname()

# A claim is a directory: mkdir either succeeds or fails, atomically, even from two servers on the shared RAID
function claim_case(name)
    isfile(joinpath(matrix_state, name * ".done")) && return false
    isfile(joinpath(matrix_state, name * ".failed")) && return false
    try
        mkdir(joinpath(matrix_state, name * ".claimed"))
    catch
        return false
    end
    write(joinpath(matrix_state, name * ".claimed", "who.txt"), "$matrix_host GPU $matrix_gpu, started $(now())\n")
    return true
end

function run_matrix_case(name)
    les_module = Module(Symbol("les_", replace(name, r"[^A-Za-z0-9]" => "_")))
    Core.eval(les_module, :(include(path) = Base.include($les_module, path)))
    empty!(ARGS); push!(ARGS, name)
    Base.include(les_module, joinpath(matrix_folder, "les.jl"))
    empty!(ARGS)
    return nothing
end

function run_matrix()
    @info "LES matrix on $matrix_host GPU $matrix_gpu: $(length(les_matrix)) cases"
    ran = String[]
    for name in les_matrix
        claim_case(name) || continue
        @info "[$matrix_host GPU $matrix_gpu] starting $name at $(now())"
        started = time()
        try
            run_matrix_case(name)
            write(joinpath(matrix_state, name * ".done"),
                  @sprintf("%s GPU %d, finished %s, %.1f h\n", matrix_host, matrix_gpu, now(), (time() - started) / 3600))
            push!(ran, name)
        catch err
            err isa InterruptException && rethrow()
            write(joinpath(matrix_state, name * ".failed"), sprint(showerror, err))
            @warn "$name failed: $(first(sprint(showerror, err), 300))"
        end
        GC.gc(true); GC.gc(true); CUDA.reclaim()
    end
    done = count(n -> isfile(joinpath(matrix_state, n * ".done")), les_matrix)
    @info "[$matrix_host GPU $matrix_gpu] nothing left to claim; this GPU ran $(length(ran)); matrix: $done of $(length(les_matrix)) done"
    return ran
end

run_matrix()
