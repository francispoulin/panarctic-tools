# Where are the hot spots? Grid + bathymetry only, on the CPU (same parameters as the script).
using NumericalEarth, Oceananigans, Oceananigans.Units, Printf
using Oceananigans.OrthogonalSphericalShellGrids: RotatedLatitudeLongitudeGrid
using Oceananigans.Grids: λnode, φnode, znode

λ₁, λ₂ = 151, 212
φ₁, φ₂ = -34.8, 27
res = 0.15
Nx, Ny, Nz = round(Int, (λ₂ - λ₁) / res), round(Int, (φ₂ - φ₁) / res), 20
depth = 4000meters
z = ExponentialDiscretization(Nz, -depth, 0; scale = depth/4, mutable = true)

grid = RotatedLatitudeLongitudeGrid(CPU(); size = (Nx, Ny, Nz), longitude = (λ₁, λ₂), latitude = (φ₁, φ₂),
                                    north_pole = (145, 0), z, halo = (7, 7, 7))

bathymetry = Metadatum(:bottom_height; dataset = ETOPO2022())
h = regrid_bathymetry(grid, bathymetry; minimum_depth = 15, interpolation_passes = 2, major_basins = 1)
hb = Array(interior(h, :, :, 1))
ℓ = (Center(), Center(), Center())

for (name, i₀, j₀) in (("fast flow", 134, 140), ("max S", 96, 359), ("min S", 1, 55))
    @printf("\n%s at (i, j) = (%d, %d) ≈ %.2f°E, %.2f°N; bottom height (m, ≥ 0 is land), i across, j down:\n",
            name, i₀, j₀, λnode(i₀, j₀, 1, grid, ℓ...), φnode(i₀, j₀, 1, grid, ℓ...))
    for j in j₀+3:-1:j₀-3
        @printf("  j = %3d ", j)
        for i in i₀-3:i₀+3
            (1 ≤ i ≤ Nx && 1 ≤ j ≤ Ny) ? @printf("%7.0f", hb[i, j]) : print("      -")
        end
        println()
    end
end

println("\nLevel depths:")
for k in Nz:-1:1
    @printf("  k = %2d   z = %7.1f m\n", k, znode(1, 1, k, grid, ℓ...))
end
