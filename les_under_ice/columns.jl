# Single-column simulations of a case in cases.jl: Oceananigans' HydrostaticFreeSurfaceModel on a (Flat, Flat,
# Bounded) grid with CATKE, the same forcing and initial state as les.jl, at 1, 4 and 8 m and on the pan-Arctic model's
# vertical grids at 50 and 100 levels (top cells ~8 and ~4 m). CATKE variants:
#   stock            CATKE as calibrated (surface TKE flux Cᵂu★ u★³ + CᵂwΔ (wΔ)³)
#   no_surface_shear Cᵂu★ = 0: no shear-driven (wave-like) surface TKE flux; the convective one is kept
#   rigid_wall       Cᵂu★ = 0 and Cˢ = 0.967 (von Kármán constant 0.4, as at the bottom): CATKE for a wall without waves
#   k_epsilon        Oceananigans' two-equation k-ε closure (TKEDissipationVerticalDiffusivity, default constants;
#                    whether they come from a calibration against LESbrary is not known here): an independent model
# Writes b (or T and S), u, v, e, κᶜ and κᵘ every 10 minutes to output/<case>_column_<variant>_<dz1|…|nz50|nz100>.jld2. Fast on a CPU:
#
#     empty!(ARGS); push!(ARGS, "free_convection")
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/columns.jl")

using Oceananigans, Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity, TKEDissipationVerticalDiffusivity
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEEquation, CATKEMixingLength
using JLD2, Printf

include(joinpath(@__DIR__, "cases.jl"))

case_name = isempty(ARGS) ? "free_convection" : ARGS[1]
case      = case_parameters(case_name)

# rigid_wall: CATKE's surface mixing-length coefficient Cˢ = 1.131 corresponds to a von Kármán constant κ = 0.47,
# larger than a rigid wall's 0.4 "because the large eddy simulations used to calibrate CATKE include surface wave
# effects" (Oceananigans, catke_mixing_length.jl); its bottom coefficient Cᵇ = 0.967 gives κ = 0.4. Under ice (a rigid
# wall, no waves): Cˢ = Cᵇ = 0.967 and no shear-driven surface TKE flux.
variants = Dict("stock"            => CATKEVerticalDiffusivity(),
                "no_surface_shear" => CATKEVerticalDiffusivity(turbulent_kinetic_energy_equation = CATKEEquation(Cᵂu★ = 0.0)),
                "rigid_wall"       => CATKEVerticalDiffusivity(mixing_length = CATKEMixingLength(Cˢ = 0.967),
                                                               turbulent_kinetic_energy_equation = CATKEEquation(Cᵂu★ = 0.0)),
                "k_epsilon"        => TKEDissipationVerticalDiffusivity())
# Uniform 1, 4, 8 m, and the pan-Arctic model's own grids at 50 and 100 levels (cases.jl, column_faces)
spacings = column_spacings

include(joinpath(@__DIR__, "column_model.jl"))   # run_column

for (variant, closure) in sort(collect(variants); by = first), spacing in spacings
    run_column(case_name, case, variant, closure, spacing)
end
