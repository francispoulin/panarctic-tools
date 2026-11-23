# NEMO's TKE vertical mixing scheme (zdftke.F90, NEMO 4.x; Gaspar et al. 1990, Blanke and Delecluse 1993, Madec et al.)
# in a standalone single column, with the same cases, forcing and output as columns.jl, so compare.jl and
# compare_lesbrary.jl show it next to CATKE. It separates what CATKE folds together (claude/catke-les-plan.md):
#   base scheme   TKE equation with shear and buoyancy production, dissipation ediss e^{3/2}/ℓ, diffusion; surface
#                 value e = ebb |τ| (also under ice, from the ice–ocean stress); mixing length √(2e)/N bounded by the
#                 distance to the boundaries (nn_mxl = 2), surface length from the stress (ln_mxl0); Prandtl number
#                 from the Richardson number (nn_pdl = 1)
#   waves         Langmuir cells (ln_lc: a TKE source over a depth set by the potential energy of the stratification
#                 against the Langmuir velocity, parameterised as 0.016 U₁₀) and penetration of surface TKE below the
#                 mixed layer (nn_etau = 1: rn_efr e_surface exp(−z/htau)); both × (1 − w(ice fraction)), with
#                 w = tanh(10 fᵢ) by default (nn_eice = 1, our "tanh")
# Defaults are NEMO's reference namelist (namzdf_tke): rn_ediff = 0.1, rn_ediss = 0.7, rn_ebb = 67.83, rn_emin = 1e-6,
# rn_emin0 = 1e-4, nn_mxl = 2, ln_mxl0 = true, nn_pdl = 1, ln_lc = true, rn_lc = 0.15, nn_etau = 1, rn_efr = 0.05,
# nn_htau = 1, nn_eice = 1; background viscosity and diffusivity 1.2e-4 and 1.2e-5 m²/s, and enhanced vertical
# diffusion of tracers where N² < 0 (namzdf: ln_zdfevd = true, nn_evdm = 0, rn_evd = 100 m²/s). NEMO convects through
# this, not through the TKE: without it (checked in a numpy port of this file, 21 Nov) free convection leaves the
# top 70 m statically unstable after 24 h; with it the convective layer is 102 m at 4 m (our LES: 107 m).
#
# Variants written by this script, for a case (CPU, seconds each):
#   nemo          the reference namelist (Langmuir cells and penetration on)
#   nemo_no_waves base scheme only (ln_lc = false, nn_etau = 0): NEMO under compact ice, or wind without waves
#   nemo_lc_only  Langmuir cells on, penetration off
# each at 1, 4, 8 m and on the pan-Arctic grids at 50 and 100 levels; buoyancy cases or seawater (T, S, TEOS-10) cases.
#
#     empty!(ARGS); push!(ARGS, "strong_wind")
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/nemo_tke.jl")
#
# Internally top-down (index 1 = surface), as NEMO: cells j = 1…N, faces k = 1…N+1 (TKE and diffusivities on faces).
# Each step, in NEMO's order: TKE (Langmuir source added, implicit solve with semi-implicit dissipation, penetration),
# mixing length and diffusivities from the new TKE, then u, v, b by implicit vertical diffusion with the surface
# fluxes, and an exact Coriolis rotation. Written from zdftke.F90 and namzdf_tke; not checked against NEMO's own 1-D
# configuration (C1D).

include(joinpath(@__DIR__, "nemo_model.jl"))   # NEMOTKEParameters, nemo_tke_column, nemo_variants

nemo_case = isempty(ARGS) ? "strong_wind" : ARGS[1]
for (label, parameters) in sort(collect(nemo_variants); by = first), spacing in column_spacings
    nemo_tke_column(nemo_case, case_parameters(nemo_case), parameters; spacing, label)
end
