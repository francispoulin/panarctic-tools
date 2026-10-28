# Transport through a strait, day by day: the model's daily 3-D output against GLORYS for the same days.
#
# The section is the row of v-faces at j through (i, j), extended east and west while there is water.
# The GLORYS transport uses GLORYS velocities (rotated into the grid's axes) on the model's own section,
# so the two differ only in how fast the water moves, not in the geometry.
#
# Usage, from the model folder (needs arctic_ocean_full.jld2 from a run):
#   julia --project=. ../panarctic-tools/panarctic-strait.jl 366 31

using NumericalEarth, Oceananigans, Oceananigans.Units, Printf, Dates
using CopernicusMarine   # GLORYS files (already downloaded by the model run)
using Oceananigans.ImmersedBoundaries: immersed_inactive_node
using Oceananigans.Operators: intrinsic_vector, Δxᶜᶠᶜ, Δzᶜᶠᶜ
using Oceananigans.Grids: λnode, φnode

i₀, j₀ = parse(Int, ARGS[1]), parse(Int, ARGS[2])

file = "arctic_ocean_full.jld2"
u = FieldTimeSeries(file, "u")
v = FieldTimeSeries(file, "v")
grid = v.grid
Nx, Ny, Nz = size(grid)

wet_face(i, k) = !immersed_inactive_node(i, j₀, k, grid, Center(), Face(), Center())
wet_column(i)  = any(k -> wet_face(i, k), 1:Nz)
wet_column(i₀) || error("no water in the v-face row through ($i₀, $j₀)")

i₁ = i₀
while i₁ > 1 && wet_column(i₁ - 1)
    i₁ -= 1
end
i₂ = i₀
while i₂ < Nx && wet_column(i₂ + 1)
    i₂ += 1
end
section = i₁:i₂

area(i, k)  = wet_face(i, k) ? Δxᶜᶠᶜ(i, j₀, k, grid) * Δzᶜᶠᶜ(i, j₀, k, grid) : 0.0
depth(i)    = sum(wet_face(i, k) ? Δzᶜᶠᶜ(i, j₀, k, grid) : 0.0 for k in 1:Nz)
width       = sum(Δxᶜᶠᶜ(i, j₀, Nz, grid) for i in section)
total_area  = sum(area(i, k) for i in section, k in 1:Nz)

ℓ = (Center(), Center(), Center())
@printf("\nSection: v-faces at j = %d, i = %d–%d (%d cells), %.2f°E %.2f°N to %.2f°E %.2f°N\n", j₀, i₁, i₂, length(section),
        λnode(i₁, j₀, 1, grid, ℓ...), φnode(i₁, j₀, 1, grid, ℓ...), λnode(i₂, j₀, 1, grid, ℓ...), φnode(i₂, j₀, 1, grid, ℓ...))
@printf("Width along the row %.1f km, cross-sectional area %.2f km², water depths [m]: %s\n",
        width / 1e3, total_area / 1e6, join([@sprintf("%.0f", depth(i)) for i in section], ", "))
println("(Transport is through the row of v-faces; positive is towards +j, i.e. towards larger j.)")

transport(vₖ) = sum(vₖ(i, k) * area(i, k) for i in section, k in 1:Nz)
surface_speed(uₖ, vₖ) = sqrt(uₖ^2 + vₖ^2)

region = BoundingBox(longitude = (-180, 180), latitude = (45, 90))
inpainting = NumericalEarth.DataWrangling.NearestNeighborInpainting(100)

println("\n  day   GLORYS mean of   model transport   GLORYS transport   model speed   GLORYS speed   (speeds at ($i₀, $j₀), surface)")
for n in eachindex(u.times)
    t = u.times[n] / days
    t < 1 && continue
    date = DateTime(1993, 1, 1) + Day(floor(Int, t - 1/2))

    uₙ, vₙ = u[n], v[n]
    model_transport = transport((i, k) -> vₙ[i, j₀, k])
    model_speed = surface_speed((uₙ[i₀, j₀, Nz] + uₙ[i₀+1, j₀, Nz]) / 2, (vₙ[i₀, j₀, Nz] + vₙ[i₀, j₀+1, Nz]) / 2)

    uᴳ = Field(Metadatum(:u_velocity; date, dataset = GLORYSDaily(), region), grid.underlying_grid; inpainting)
    vᴳ = Field(Metadatum(:v_velocity; date, dataset = GLORYSDaily(), region), grid.underlying_grid; inpainting)
    vᴳ_grid(i, k) = intrinsic_vector(i, j₀, k, grid, uᴳ[i, j₀, k], vᴳ[i, j₀, k])[2]
    glorys_transport = transport(vᴳ_grid)
    glorys_speed = surface_speed(uᴳ[i₀, j₀, Nz], vᴳ[i₀, j₀, Nz])

    @printf("  %3.0f   %-14s   %+8.3f Sv       %+8.3f Sv          %5.2f m/s     %5.2f m/s\n", t,
            Dates.format(date, "d u yyyy"), model_transport / 1e6, glorys_transport / 1e6, model_speed, glorys_speed)
end
