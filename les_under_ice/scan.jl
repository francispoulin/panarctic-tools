# Step 2 of claude/catke-les-plan.md: shear alone. Which CATKE parameters carry the wave effect? Scans the surface
# mixing-length coefficient Cˢ and the shear-driven surface TKE flux coefficient Cᵂu★ in single columns, against our
# LES of a case (default strong_wind: wind without waves, the situation under ice). Errors at 6, 12, 24 h, between 4 and
# 160 m (as the calibration), the LES averaged onto the column's cells:
#   b error  ‖b − b_LES‖ / ‖b_LES − b_initial‖        u error  ‖(u, v) − (u, v)_LES‖ / ‖(u, v)_LES‖
# Writes output/scan_<case>.md (every pair, sorted by b error) and output/scan_<case>.png (b-error map). On a CPU:
#
#     empty!(ARGS); push!(ARGS, "strong_wind")                 # second argument: Δz (default 4)
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/scan.jl")
#
# Stock CATKE: Cˢ = 1.131 (von Kármán constant 0.47, wave-influenced), Cᵂu★ = 3.179. Rigid wall: Cˢ = 0.967 (κ = 0.4).
# About 50 columns of a few seconds each.

using Oceananigans, Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEEquation, CATKEMixingLength
using Oceananigans.OutputReaders: FieldTimeSeries
using CairoMakie, JLD2, Printf, Statistics

include(joinpath(@__DIR__, "cases.jl"))

case_name = isempty(ARGS) ? "strong_wind" : ARGS[1]
Δz_scan   = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 4.0
case      = case_parameters(case_name)

Cˢ_values   = [0.6, 0.7, 0.8, 0.9, 0.967, 1.05, 1.131, 1.25]
Cᵂu★_values = [0.0, 0.5, 1.0, 2.0, 3.179]
check_times = [6, 12, 24] .* Oceananigans.Units.hour

#####
##### The LES, averaged onto the column grid
#####

les_file = joinpath(output_folder(), case_name * "_les.jld2")
isfile(les_file) || error("No LES $les_file: run les.jl $case_name first")

function les_profiles(file, zᶠ)
    B, U, V = FieldTimeSeries(file, "b"), FieldTimeSeries(file, "u"), FieldTimeSeries(file, "v")
    z = collect(znodes(B.grid, Center()))
    average(values) = [mean(values[(z .> zᶠ[k]) .& (z .< zᶠ[k+1])]) for k in 1:length(zᶠ)-1]
    column(fts, n) = average(vec(Array(interior(fts[n], 1, 1, :))))
    index(t) = argmin(abs.(B.times .- t))
    return Dict(t => (b = column(B, index(t)), u = column(U, index(t)), v = column(V, index(t))) for t in check_times)
end

column_grid() = RectilinearGrid(CPU(); size = round(Int, domain_depth / Δz_scan), z = (-domain_depth, 0),
                                topology = (Flat, Flat, Bounded))
grid₀  = column_grid()
zᶠ     = collect(znodes(grid₀, Face()))
zc     = collect(znodes(grid₀, Center()))
les    = les_profiles(les_file, zᶠ)
window = (zc .> -160) .& (zc .< -4)
b₀     = three_layer_buoyancy.(zc)

#####
##### One column: errors at the check times
#####

function column_errors(closure)
    grid = column_grid()
    b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(case.Jb),
                                    bottom = GradientBoundaryCondition(three_layer.N²_deep))
    u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(case.τx))
    model = HydrostaticFreeSurfaceModel(grid; closure, coriolis = FPlane(f = case.f), tracers = :b,
                                        buoyancy = BuoyancyTracer(), boundary_conditions = (u = u_bcs, b = b_bcs))
    set!(model, b = three_layer_buoyancy, e = 1e-6)
    simulation = Simulation(model; Δt = 1minute, stop_time = 0.0, verbose = false)
    b_errors, u_errors = Float64[], Float64[]
    for t in check_times
        simulation.stop_time = t
        run!(simulation)
        b = vec(interior(model.tracers.b, 1, 1, :))
        u = vec(interior(model.velocities.u, 1, 1, :))
        v = vec(interior(model.velocities.v, 1, 1, :))
        ref = les[t]
        push!(b_errors, sqrt(sum(abs2, (b .- ref.b)[window]) / sum(abs2, (ref.b .- b₀)[window])))
        push!(u_errors, sqrt(sum(abs2, (u .- ref.u)[window]) + sum(abs2, (v .- ref.v)[window])) /
                        sqrt(sum(abs2, ref.u[window]) + sum(abs2, ref.v[window])))
    end
    return (b = mean(b_errors), u = mean(u_errors), b_end = b_errors[end])
end

closure_for(Cˢ, Cᵂu★) = CATKEVerticalDiffusivity(mixing_length = CATKEMixingLength(; Cˢ),
                                                 turbulent_kinetic_energy_equation = CATKEEquation(; Cᵂu★))

results = NamedTuple[]
for Cˢ in Cˢ_values, Cᵂu★ in Cᵂu★_values
    errors = column_errors(closure_for(Cˢ, Cᵂu★))
    push!(results, (; Cˢ, Cᵂu★, errors...))
    @info @sprintf("Cˢ = %.3f, Cᵂu★ = %.3f: b error %.1f %% (24 h: %.1f %%), u error %.1f %%",
                   Cˢ, Cᵂu★, 100errors.b, 100errors.b_end, 100errors.u)
end

#####
##### Table and figure
#####

sorted = sort(results; by = r -> r.b)
stock  = only(filter(r -> r.Cˢ == 1.131 && r.Cᵂu★ == 3.179, results))
io = IOBuffer()
println(io, "# Scan of Cˢ and Cᵂu★ against our LES of $case_name (Δz = $(Δz_scan) m)\n")
println(io, "Errors averaged over 6, 12, 24 h, between 4 and 160 m depth. Stock CATKE: Cˢ = 1.131, Cᵂu★ = 3.179: ",
            @sprintf("b error %.1f %%, u error %.1f %%.", 100stock.b, 100stock.u), " Rigid wall: Cˢ = 0.967 (κ = 0.4).\n")
println(io, "| Cˢ | Cᵂu★ | b error | b error at 24 h | u error |")
println(io, "|---|---|---|---|---|")
for r in sorted
    println(io, @sprintf("| %.3f | %.3f | %.1f %% | %.1f %% | %.1f %% |", r.Cˢ, r.Cᵂu★, 100r.b, 100r.b_end, 100r.u))
end
text = String(take!(io))
print(text)
write(joinpath(output_folder(), "scan_$(case_name).md"), text)

fig = Figure(size = (900, 600))
ax = Axis(fig[1, 1]; xlabel = "Cˢ (surface mixing-length coefficient)", ylabel = "Cᵂu★ (surface shear TKE flux)",
          title = "$case_name: b error of CATKE columns (Δz = $(Δz_scan) m) against our LES")
errors = [only(filter(r -> r.Cˢ == Cˢ && r.Cᵂu★ == Cᵂu★, results)).b for Cˢ in Cˢ_values, Cᵂu★ in Cᵂu★_values]
hm = heatmap!(ax, 1:length(Cˢ_values), 1:length(Cᵂu★_values), 100 .* errors; colormap = :viridis)
ax.xticks = (1:length(Cˢ_values), string.(Cˢ_values))
ax.yticks = (1:length(Cᵂu★_values), string.(Cᵂu★_values))
for i in eachindex(Cˢ_values), j in eachindex(Cᵂu★_values)
    text!(ax, i, j; text = @sprintf("%.0f", 100errors[i, j]), align = (:center, :center), color = :white)
end
Colorbar(fig[1, 2], hm; label = "b error (%)")
save(joinpath(output_folder(), "scan_$(case_name).png"), fig)
@info "Wrote output/scan_$(case_name).md and .png"
