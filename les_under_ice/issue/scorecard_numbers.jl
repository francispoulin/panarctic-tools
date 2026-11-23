# Numbers for the Oceananigans issue's scorecards, recomputed from the column and LES files already in output/, all with
# the current depth measure (cases.jl case_depth, interpolated between levels since 23 Nov). The 22 Nov scorecard was
# typed in from tables made before that change; the ice cases (10–26 m deep, 4 m cells) moved most.
#
# Depth error (h_column − h_LES) / h_LES and the depths themselves, for:
#   scorecard 1 (three schemes, 4 m, Δt = 1 min):  CATKE stock, k-ε, NEMO (wave terms off without waves, on with waves)
#                                                  from output/<case>_column_<variant>_dz4.jld2 (columns.jl, nemo_model.jl)
#   scorecard 2 (adapting CATKE, Δt = 10 min):     stock and cap0.1_kappa041_flux_half on dz4 and glorys50
#                                                  from output/ice_variants/ (ice_variants.jl)
# A missing file gives NaN. Writes output/scorecard_numbers.csv and prints it. Reads files only (a minute or two):
#
#     include("/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/les_under_ice/issue/scorecard_numbers.jl")

using Oceananigans, Oceananigans.Units
using JLD2, Printf, Statistics

include(joinpath(@__DIR__, "..", "cases.jl"))

function scorecard_coarse_grain(values, z_fine, zᶠ)
    return map(1:length(zᶠ)-1) do k
        inside = (z_fine .> zᶠ[k]) .& (z_fine .< zᶠ[k+1])
        any(inside) ? mean(values[inside]) : values[argmin(abs.(z_fine .- (zᶠ[k] + zᶠ[k+1]) / 2))]
    end
end

# (case, has waves)
scorecard_cases = [
    (convection_case_name(1e-5, 1.2e-7), false), (convection_case_name(1e-5, 2.4e-7), false),
    (convection_case_name(1e-5, 4.8e-7), false), (convection_case_name(4e-5, 4.8e-7), false),
    (convection_case_name(2.5e-6, 1.2e-7), false), ("arctic_brine_convection", false),
    ("kato_phillips", false), ("rotating_stress", false), ("arctic_ice_drag", false), ("arctic_ice_drag_brine", false),
    ("arctic_ice_drag_melt", false),
    ("weak_wind_strong_cooling", false), ("med_wind_med_cooling", false), ("strong_wind_weak_cooling", false),
    ("strong_wind", false),
    ("weak_wind_strong_cooling_stokes", true), ("med_wind_med_cooling_stokes", true),
    ("strong_wind_weak_cooling_stokes", true), ("strong_wind_stokes", true)]

# (column label, file) for a case
function scorecard_files(name, waves, folder)
    nemo = waves ? "nemo" : "nemo_no_waves"
    nemo_file = joinpath(folder, "$(name)_column_$(nemo)_dz4.jld2")
    waves || isfile(nemo_file) || (nemo_file = joinpath(folder, "$(name)_column_nemo_dz4.jld2"))   # convection: same
    iv(variant, grid) = joinpath(folder, "ice_variants", "$(name)_$(variant)_$(grid)_dt10.jld2")
    return [("catke_dz4_dt1",               joinpath(folder, "$(name)_column_stock_dz4.jld2")),
            ("k_epsilon_dz4_dt1",           joinpath(folder, "$(name)_column_k_epsilon_dz4.jld2")),
            ("nemo_dz4_dt1",                nemo_file),
            ("catke_stock_dz4_dt10",        iv("stock", "dz4")),
            ("catke_adapted_dz4_dt10",      iv("cap0.1_kappa041_flux_half", "dz4")),
            ("catke_stock_glorys50_dt10",   iv("stock", "glorys50")),
            ("catke_adapted_glorys50_dt10", iv("cap0.1_kappa041_flux_half", "glorys50"))]
end

function scorecard_numbers(folder)
    lines = ["case,column,h_les_m,h_column_m,error"]
    for (name, waves) in scorecard_cases
        les_file = joinpath(folder, name * "_les.jld2")
        isfile(les_file) || (@warn "No LES for $name"; continue)
        les = les_profiles(les_file)
        case = case_parameters(name)
        for (label, file) in scorecard_files(name, waves, folder)
            if !isfile(file)
                push!(lines, "$name,$label,NaN,NaN,NaN")
                continue
            end
            col = column_profiles(file)
            b_les = scorecard_coarse_grain(les.b[argmin(abs.(les.t .- col.t[end]))], les.zc, col.zf)
            h_les = case_depth(name, case, b_les, col.zc)
            h_col = case_depth(name, case, col.b[end], col.zc)
            push!(lines, @sprintf("%s,%s,%.2f,%.2f,%.4f", name, label, h_les, h_col, (h_col - h_les) / h_les))
        end
    end
    text = join(lines, "\n") * "\n"
    write(joinpath(folder, "scorecard_numbers.csv"), text)
    print(text)
    return text
end

scorecard_numbers(output_folder())
@info "Wrote $(joinpath(output_folder(), "scorecard_numbers.csv"))"
