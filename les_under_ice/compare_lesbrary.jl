# Step 1 of claude/catke-les-plan.md: validate our LES and the CATKE columns against the calibration LES of Wagner et
# al. (LESbrary profiles from github.com/glwagner/SingleColumnModelCalibration.jl, downloaded to data/profiles/).
#
# For a case of their 24-hour suite, plots b and u at 6, 12 and 24 h from: their LES at 1, 2 and 4 m; our LES
# (les.jl; their strong_wind has Stokes drift, so it is compared with our strong_wind_stokes); and our CATKE columns
# (columns.jl). Writes a table of boundary-layer depths and b errors against their 1 m LES to
# output/lesbrary_<case>.md and the figure to output/lesbrary_<case>.png. In a REPL (CPU is enough):
#
#     empty!(ARGS); push!(ARGS, "strong_wind")                 # or "free_convection"
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/compare_lesbrary.jl")
#
# Their files are read from the raw JLD2 contents (timeseries/<field>/<iteration>, grid/...), not through
# Oceananigans' FieldTimeSeries: they were written by an older Oceananigans whose serialized grid may not load.

using Oceananigans, Oceananigans.Units
using Oceananigans.OutputReaders: FieldTimeSeries
using CairoMakie, JLD2, Printf, Statistics

include(joinpath(@__DIR__, "cases.jl"))

their_case = isempty(ARGS) ? "strong_wind" : ARGS[1]
hours      = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 24
our_les_case = Dict("strong_wind" => "strong_wind_stokes", "free_convection" => "free_convection")
our_case   = get(our_les_case, their_case, their_case)
folder     = output_folder()
profiles   = joinpath(@__DIR__, "data", "profiles", "$(hours)_hour_suite")

#####
##### Their LES, from the raw JLD2
#####

interior_values(a) = vec(collect(parent(a)))

# Strip halo points if the profile carries them
function strip_halos(values, Nz, Hz)
    n = length(values)
    n == Nz          && return values
    n == Nz + 2Hz    && return values[Hz+1:Hz+Nz]
    n == Nz + 2Hz + 1 && return values[Hz+1:Hz+Nz+1]   # a face field with halos
    error("Profile of length $n with Nz = $Nz and Hz = $Hz")
end

function read_lesbrary(file)
    jldopen(file) do f
        Nz, Hz = f["grid/Nz"], f["grid/Hz"]
        grid_keys = keys(f["grid"])
        zkey = findfirst(k -> k in grid_keys, ["zᵃᵃᶜ", "zC", "zc"])
        zkey === nothing && error("No z-centres in $file; grid keys: $grid_keys")
        zc = strip_halos(interior_values(f["grid/" * ["zᵃᵃᶜ", "zC", "zc"][zkey]]), Nz, Hz)
        iterations = sort(filter(!isnothing, tryparse.(Int, keys(f["timeseries/t"]))))   # skips "serialized" entries
        t = [f["timeseries/t/$i"] for i in iterations]
        field(name) = [strip_halos(interior_values(f["timeseries/$name/$i"]), Nz, Hz) for i in iterations]
        parameters = Dict(k => f["parameters/$k"] for k in keys(f["parameters"]))
        return (; zc, t, b = field("b"), u = field("u"), parameters)
    end
end

function their_file(resolution)
    dir = joinpath(profiles, resolution)
    isdir(dir) || return nothing
    matches = filter(n -> startswith(n, their_case * "_") && endswith(n, "statistics.jld2"), readdir(dir))
    preferred = filter(n -> occursin("with_tracer_instantaneous", n), matches)
    name = isempty(preferred) ? (isempty(matches) ? nothing : first(matches)) : first(preferred)
    return name === nothing ? nothing : joinpath(dir, name)
end

theirs = Dict{String, Any}()
for resolution in ("1m", "2m", "4m")
    file = their_file(resolution)
    file === nothing && continue
    theirs[resolution] = read_lesbrary(file)
end
haskey(theirs, "1m") || error("No $their_case profile at 1m in $profiles")
reference = theirs["1m"]

println("Their $their_case ($(hours) h suite) parameters:")
for k in sort(collect(keys(reference.parameters)))
    v = reference.parameters[k]
    v isa Union{Number, AbstractString} && println("  $k = $v")
end
case = haskey(cases, our_case) ? cases[our_case] : nothing
case === nothing || println("Our $our_case: τx = $(case.τx), Jb = $(case.Jb), f = $(case.f), Stokes drift: $(case.stokes)")

#####
##### Ours
#####

profile(fts, n) = vec(Array(interior(fts[n], 1, 1, :)))
nearest(times, t) = argmin(abs.(times .- t))

ours = Dict{String, Any}()
les_file = joinpath(folder, our_case * "_les.jld2")
if isfile(les_file)
    B, U = FieldTimeSeries(les_file, "b"), FieldTimeSeries(les_file, "u")
    ours["our LES (2 m)"] = (; zc = collect(znodes(B.grid, Center())), t = B.times,
                              b = [profile(B, n) for n in eachindex(B.times)], u = [profile(U, n) for n in eachindex(U.times)])
end
for name in sort(filter(n -> startswith(n, our_case * "_column_") && endswith(n, ".jld2"), readdir(folder)))
    B, U = FieldTimeSeries(joinpath(folder, name), "b"), FieldTimeSeries(joinpath(folder, name), "u")
    label = "CATKE " * replace(name, our_case * "_column_" => "", ".jld2" => "")
    ours[label] = (; zc = collect(znodes(B.grid, Center())), t = B.times,
                     b = [profile(B, n) for n in eachindex(B.times)], u = [profile(U, n) for n in eachindex(U.times)])
end

#####
##### Depths and errors against their 1 m LES, at the end
#####

function boundary_layer_depth(b, z; Δb = 1e-5, z_ref = -10.0)
    reference_b = b[argmin(abs.(z .- z_ref))]
    k = findlast(k -> z[k] < z_ref && b[k] < reference_b - Δb, eachindex(b))
    return k === nothing ? -z[1] : -z[k]
end

# Their 1 m profile interpolated onto other centres (linear)
function interpolate(values, z_from, z_to)
    order = sortperm(z_from)
    zs, vs = z_from[order], values[order]
    return map(z_to) do z
        z ≤ zs[1] && return vs[1]
        z ≥ zs[end] && return vs[end]
        k = searchsortedlast(zs, z)
        w = (z - zs[k]) / (zs[k+1] - zs[k])
        (1 - w) * vs[k] + w * vs[k+1]
    end
end

t_end = min(hours * hour, reference.t[end])
n_ref = nearest(reference.t, t_end)
b_ref_end = reference.b[n_ref]
b_ref_start = reference.b[1]
h_ref = boundary_layer_depth(b_ref_end, reference.zc)

rows = String[]
entries = vcat([("LESbrary " * r, theirs[r]) for r in ("1m", "2m", "4m") if haskey(theirs, r)],
               [(k, ours[k]) for k in sort(collect(keys(ours)))])
for (label, data) in entries
    n = nearest(data.t, t_end)
    b = data.b[n]
    h = boundary_layer_depth(b, data.zc)
    on_grid = interpolate(b_ref_end, reference.zc, data.zc)
    change  = on_grid .- interpolate(b_ref_start, reference.zc, data.zc)
    top = (data.zc .> -160) .& (data.zc .< -4)   # as the calibration: without the top 4 m
    err = sqrt(sum(abs2, (b .- on_grid)[top]) / sum(abs2, change[top]))
    push!(rows, @sprintf("| %s | %.0f %% | %.1f | %+.1f |", label, 100err, h, h - h_ref))
end

io = IOBuffer()
println(io, "# $their_case ($(hours)-hour suite): our LES and CATKE columns against the calibration LES\n")
println(io, "Reference: LESbrary 1 m. At $(round(Int, t_end / hour)) h. b error: ‖b − b_ref‖ / ‖b_ref − b_ref(0)‖ between 4 and ",
            "160 m depth, on each profile's own levels. h: first depth below 10 m where b is below its 10 m value by ",
            "more than 1e-5 m s⁻². Our LES/columns for this case: $our_case.\n")
println(io, "| Profile | b error | h (m) | h − h(LESbrary 1 m) (m) |")
println(io, "|---|---|---|---|")
foreach(r -> println(io, r), rows)
text = String(take!(io))
print(text)
write(joinpath(folder, "lesbrary_$(their_case).md"), text)

#####
##### Figure
#####

times = [t for t in (6hours, 12hours, 24hours) if t ≤ reference.t[end] + 1]
fig = Figure(size = (1500, 900))
colors = Makie.wong_colors()
for (i, t) in enumerate(times)
    axb = Axis(fig[1, i]; title = "b at $(round(Int, t / hour)) h", xlabel = "b (m s⁻²)", ylabel = i == 1 ? "z (m)" : "")
    axu = Axis(fig[2, i]; title = "u at $(round(Int, t / hour)) h", xlabel = "u (m/s)", ylabel = i == 1 ? "z (m)" : "")
    for (c, (label, data)) in enumerate(entries)
        n = nearest(data.t, t)
        style = startswith(label, "LESbrary") ? :solid : startswith(label, "our LES") ? :dash : :dot
        width = startswith(label, "LESbrary 1m") ? 3 : 1.5
        lines!(axb, data.b[n], data.zc; color = colors[mod1(c, 7)], linestyle = style, linewidth = width, label)
        lines!(axu, data.u[n], data.zc; color = colors[mod1(c, 7)], linestyle = style, linewidth = width, label)
    end
    ylims!(axb, -200, 0); ylims!(axu, -200, 0)
    i == 1 && axislegend(axb; position = :rb, labelsize = 9)
end
Label(fig[0, :], "$their_case ($(hours)-hour suite): LESbrary (solid), our LES (dashed), CATKE columns (dotted)"; fontsize = 18)
save(joinpath(folder, "lesbrary_$(their_case).png"), fig)
@info "Wrote output/lesbrary_$(their_case).md and .png"
