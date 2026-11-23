# Step 3 of the CATKE recalibration plan (claude/catke-les-plan.md): which of CATKE's 24 active parameters matter?
#
# Morris screening ("elementary effects"; Morris 1991, Campolongo et al. 2007): from r random starting points in the
# parameter space, change the parameters one at a time by a fixed step and record the change of the error. For each
# parameter: μ* = mean |effect| (how much it matters), σ = spread of the effect (whether it acts nonlinearly or through
# other parameters). r (24 + 1) parameter sets × the cases below, single columns at 4 m, CPU, spread over worker
# processes (Distributed, one thread each).
#
# Parameters: the 24 non-zero constants of CATKEMixingLength and CATKEEquation (Cᵇ, the bottom wall, is left out), each
# varied log-uniformly between a third of and three times its default. Morris grid of 4 levels, step 2/3 of the range.
#
# Error of a run: the mixed-layer depth error |h − h_LES| / h_LES (depth interpolated between levels, cases.jl
# case_depth), against the LES coarse-grained to the column. Groups: targets (no waves, CATKE too deep) and guards
# (CATKE right). A parameter worth calibrating has a large μ* on the targets; one that also has a large μ* on the
# guards must be calibrated with the guards in the loss.
#
# Writes output/screen.md (ranking), output/screen.png (μ* targets vs guards) and output/screen.jld2 (every parameter
# set and every error, for step 4). In a fresh REPL (it starts the workers itself; ~30–60 min):
#
#     empty!(ARGS); push!(ARGS, "24", "10")         # workers, trajectories r (defaults 24 and 10)
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/screen.jl")

using Distributed

screen_workers      = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 24
screen_trajectories = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 10
empty!(ARGS)

screen_folder = @__DIR__
existing_workers = workers() == [1] ? 0 : nworkers()
existing_workers < screen_workers &&
    addprocs(screen_workers - existing_workers; exeflags = "--project=$(dirname(Base.active_project()))",
             env = ["JULIA_NUM_THREADS" => "1"])

@everywhere begin
    using Logging
    using Oceananigans, Oceananigans.Units
    using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
    using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEEquation, CATKEMixingLength
    include(joinpath($screen_folder, "cases.jl"))
    include(joinpath($screen_folder, "column_model.jl"))
end

using CairoMakie, JLD2, Printf, Random, Statistics

# The parameters screened: (name, default, which constructor)
screen_parameters = [
    (:Cˢ, 1.131, :ml), (:Cʰⁱu, 0.242, :ml), (:Cʰⁱc, 0.098, :ml), (:Cʰⁱe, 0.548, :ml), (:Cˡᵒu, 0.361, :ml),
    (:Cˡᵒc, 0.369, :ml), (:Cˡᵒe, 7.863, :ml), (:Cᵘⁿu, 0.37, :ml), (:Cᵘⁿc, 0.572, :ml), (:Cᵘⁿe, 1.447, :ml),
    (:Cᶜu, 3.705, :ml), (:Cᶜc, 4.793, :ml), (:Cᶜe, 3.642, :ml), (:Cᵉc, 0.112, :ml), (:Cˢᵖ, 0.505, :ml),
    (:CRiᵟ, 1.02, :ml), (:CRi⁰, 0.254, :ml),
    (:CʰⁱD, 0.579, :eq), (:CˡᵒD, 1.604, :eq), (:CᵘⁿD, 0.923, :eq), (:CᶜD, 3.254, :eq), (:Cᵂu★, 3.179, :eq),
    (:CᵂwΔ, 0.383, :eq), (:Cᵂϵ, 1.0, :eq)]
k_parameters = length(screen_parameters)

screen_targets = ["rotating_stress", "arctic_ice_drag", "arctic_ice_drag_melt", "med_wind_med_cooling", "strong_wind"]
screen_guards  = ["kato_phillips", convection_case_name(1e-5, 2.4e-7), "arctic_brine_convection",
                  "med_wind_med_cooling_stokes"]
screen_cases   = vcat(screen_targets, screen_guards)
screen_spacing = 4

# Unit coordinate x ∈ [0, 1] → parameter value, log-uniform between default/3 and 3 default
parameter_value(x, default) = default * 3.0^(2x - 1)

# Morris trajectories on a 4-level grid {0, 1/3, 2/3, 1} with step Δ = 2/3: each is k + 1 points, the parameters
# changed one at a time in random order, each by +Δ or −Δ (whichever stays in [0, 1])
function morris_trajectories(k, r; rng = MersenneTwister(2026))
    Δ = 2 / 3
    trajectories = Vector{Tuple{Matrix{Float64}, Vector{Int}}}()
    for _ in 1:r
        x = rand(rng, [0, 1 / 3], k)                     # base points from which +Δ stays in range
        order = randperm(rng, k)
        points = zeros(k + 1, k)
        points[1, :] = x
        for (step, i) in enumerate(order)
            x = copy(x)
            x[i] += Δ
            points[step + 1, :] = x
        end
        push!(trajectories, (points, order))
    end
    return trajectories, Δ
end

@everywhere function screened_closure(names, values, constructors)
    ml = NamedTuple{Tuple(n for (n, c) in zip(names, constructors) if c == :ml)}(Tuple(v for (v, c) in zip(values, constructors) if c == :ml))
    eq = NamedTuple{Tuple(n for (n, c) in zip(names, constructors) if c == :eq)}(Tuple(v for (v, c) in zip(values, constructors) if c == :eq))
    return CATKEVerticalDiffusivity(mixing_length = CATKEMixingLength(; ml...),
                                    turbulent_kinetic_energy_equation = CATKEEquation(; eq...))
end

# Returns the final buoyancy profile, or the error text (a parameter set that blows up is recorded as missing, but the
# error is kept so that a bug that fails every run is reported instead of becoming NaN everywhere)
@everywhere function screened_run(case_name, names, values, constructors, spacing)
    case = case_parameters(case_name)
    b = with_logger(NullLogger()) do
        try
            final_buoyancy(case, screened_closure(names, values, constructors), spacing)
        catch err
            sprint(showerror, err)
        end
    end
    return b
end

function screen_coarse_grain(values, z_fine, zᶠ)
    return map(1:length(zᶠ)-1) do k
        inside = (z_fine .> zᶠ[k]) .& (z_fine .< zᶠ[k+1])
        any(inside) ? mean(values[inside]) : values[argmin(abs.(z_fine .- (zᶠ[k] + zᶠ[k+1]) / 2))]
    end
end

function run_screen(folder, r)
    trajectories, Δ = morris_trajectories(k_parameters, r)
    names        = first.(screen_parameters)
    defaults     = [p[2] for p in screen_parameters]
    constructors = last.(screen_parameters)
    points = reduce(vcat, first.(trajectories))                   # all parameter sets, (r (k+1)) × k
    values_of(x) = parameter_value.(x, defaults)

    # The LES depth of each case on the 4 m column grid
    zᶠ = Dict(name => column_faces(screen_spacing, case_parameters(name).depth) for name in screen_cases)
    zc = Dict(name => (zᶠ[name][1:end-1] .+ zᶠ[name][2:end]) ./ 2 for name in screen_cases)
    h_les = Dict(map(screen_cases) do name
        les = les_profiles(joinpath(folder, name * "_les.jld2"))
        b = screen_coarse_grain(les.b[end], les.zc, zᶠ[name])
        name => case_depth(name, case_parameters(name), b, zc[name])
    end)

    # Everything a worker needs is computed here and sent as plain data: the function given to pmap must not call
    # functions or read globals that exist only on this process (Julia ≥ 1.12 refuses to rebuild it on the worker:
    # "UndefVarError: `#parameter_value` not defined in `Main`", 23 Nov)
    spacing = screen_spacing
    parameter_sets = [values_of(points[i, :]) for i in 1:size(points, 1)]
    jobs = [(name, parameter_sets[i]) for name in screen_cases for i in 1:size(points, 1)]
    job_keys = [(name, i) for name in screen_cases for i in 1:size(points, 1)]

    # Fail fast: one column with the default parameters on a worker, errors not caught (~1 min)
    check = remotecall_fetch(screened_run, first(workers()), "kato_phillips", names, defaults, constructors, spacing)
    check isa AbstractString && error("Default-parameter column failed on a worker:\n$check")
    @info "Check passed (default parameters, kato_phillips, on worker $(first(workers())))"

    @info "Screening: $(length(jobs)) columns on $(nworkers()) workers"
    profiles = pmap(jobs) do (name, parameters)
        screened_run(name, names, parameters, constructors, spacing)
    end
    failures = [b for b in profiles if b isa AbstractString]
    if !isempty(failures)
        @warn "$(length(failures)) of $(length(jobs)) columns failed; distinct errors:"
        foreach(m -> println("  ", first(m, 300)), unique(failures))
    end
    errors = Dict{Tuple{String, Int}, Float64}()
    for ((name, i), b) in zip(job_keys, profiles)
        errors[(name, i)] = b isa AbstractString ? NaN :
                            abs(case_depth(name, case_parameters(name), b, zc[name]) - h_les[name]) / h_les[name]
    end

    # Elementary effects per parameter and case, then μ* and σ over trajectories and cases of each group
    effects = Dict{Tuple{String, Int}, Vector{Float64}}()          # (case, parameter) => effects
    offset = 0
    for (traj, order) in trajectories
        for name in screen_cases, (step, i) in enumerate(order)
            e = (errors[(name, offset + step + 1)] - errors[(name, offset + step)]) / Δ
            push!(get!(effects, (name, i), Float64[]), e)
        end
        offset += k_parameters + 1
    end
    group_stats(group, i) = (all_effects = filter(isfinite, reduce(vcat, [effects[(c, i)] for c in group]));
                             (mean(abs.(all_effects)), std(all_effects), mean(all_effects)))
    ranking = sort(1:k_parameters; by = i -> -group_stats(screen_targets, i)[1])

    io = IOBuffer()
    println(io, "# CATKE parameter screening (Morris elementary effects)\n")
    println(io, "$(r) trajectories × $(k_parameters + 1) parameter sets × $(length(screen_cases)) cases at $(screen_spacing) m. ",
                "Each parameter varied log-uniformly in [default/3, 3 default]; error = |h − h_LES| / h_LES. ",
                "μ*: mean |effect| (importance); σ: spread (nonlinear or interacting); mean: sign of the effect ",
                "(> 0: raising the parameter increases the error). Blown-up runs: $(count(isnan, values(errors))).\n")
    println(io, "Targets: $(join(screen_targets, ", ")). Guards: $(join(screen_guards, ", ")).\n")
    println(io, "| Parameter | default | μ* targets | σ targets | mean targets | μ* guards | σ guards |")
    println(io, "|---|---|---|---|---|---|---|")
    for i in ranking
        μt, σt, mt = group_stats(screen_targets, i)
        μg, σg, _  = group_stats(screen_guards, i)
        println(io, @sprintf("| %s | %.3f | %.3f | %.3f | %+.3f | %.3f | %.3f |", names[i], defaults[i], μt, σt, mt, μg, σg))
    end
    println(io, "\nμ* on the targets, per case:\n")
    println(io, "| Parameter | " * join(screen_targets, " | ") * " |")
    println(io, "|---|" * repeat("---|", length(screen_targets)))
    for i in ranking[1:min(10, end)]
        cells = [@sprintf("%.3f", mean(abs.(filter(isfinite, effects[(c, i)])))) for c in screen_targets]
        println(io, "| $(names[i]) | " * join(cells, " | ") * " |")
    end
    text = String(take!(io))
    print(text)
    write(joinpath(folder, "screen.md"), text)

    fig = Figure(size = (900, 750))
    ax = Axis(fig[1, 1]; xlabel = "μ* on the targets (no waves)", ylabel = "μ* on the guards",
              title = "CATKE parameters: importance under ice against importance where CATKE is already right")
    μt = [group_stats(screen_targets, i)[1] for i in 1:k_parameters]
    μg = [group_stats(screen_guards, i)[1] for i in 1:k_parameters]
    scatter!(ax, μt, μg; markersize = 10)
    text!(ax, μt, μg; text = String.(names), offset = (5, 3), fontsize = 13)
    save(joinpath(folder, "screen.png"), fig)

    jldsave(joinpath(folder, "screen.jld2"); names = String.(names), defaults, points = values_of.(eachrow(points)),
            cases = screen_cases, errors = [errors[(c, i)] for c in screen_cases, i in 1:size(points, 1)],
            h_les = [h_les[c] for c in screen_cases])
    return text
end

run_screen(output_folder(), screen_trajectories)
@info "Wrote $(joinpath(output_folder(), "screen.md")), screen.png and screen.jld2"
