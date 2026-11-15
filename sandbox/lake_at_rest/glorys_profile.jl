# One GLORYS temperature and salinity profile for the lake-at-rest test: the deepest water column within 0.5° of
# 0°E 72°N (Mohn Ridge, where the 2005 run had a 3 m/s jet at 2000 m) on 1 September 2005.
# Writes glorys_profile_0E_72N.csv next to this file: depth (m, positive down), T (°C), S (g/kg).
#
# Run with the nested environment of PanArcticShenanigans (it has NumericalEarth and the GLORYS download):
#   julia --project=<PanArcticShenanigans>/nested glorys_profile.jl

using NumericalEarth, CopernicusMarine, Dates, Printf
using NumericalEarth.DataWrangling: BoundingBox, Metadatum, default_download_directory, NCDatasets, Downloads
using NumericalEarth.DataWrangling.GLORYS: GLORYSDaily, GLORYS_dataset_variable_names
using .NCDatasets: Dataset

region = BoundingBox(longitude = (-0.5, 0.5), latitude = (71.5, 72.5))
date   = DateTime(2005, 9, 1)
dir    = default_download_directory(GLORYSDaily())

function read_column_variable(name)
    path = Downloads.download(Metadatum(name; dataset = GLORYSDaily(), date, region, dir))
    ds = Dataset(path)
    data  = Float64.(coalesce.(ds[GLORYS_dataset_variable_names[name]][:, :, :, 1], NaN))
    depth = Float64.(ds["depth"][:])
    close(ds)
    return data, depth
end

T, depth = read_column_variable(:temperature)
S, _     = read_column_variable(:salinity)

wet_levels = [count(isfinite, T[i, j, :]) for i in axes(T, 1), j in axes(T, 2)]
i, j = Tuple(argmax(wet_levels))
Nwet = wet_levels[i, j]
@info @sprintf("Deepest column: %d wet levels, down to %.0f m", Nwet, depth[Nwet])

file = joinpath(@__DIR__, "glorys_profile_0E_72N.csv")
open(file, "w") do io
    println(io, "depth_m,temperature_C,salinity_gkg")
    for k in 1:Nwet
        @printf(io, "%.2f,%.4f,%.4f\n", depth[k], T[i, j, k], S[i, j, k])
    end
end
@info "Wrote $file"
