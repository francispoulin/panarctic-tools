# Gravitational stability of the initial state (t = 0) of the PanArcticShenanigans run.
#
# Builds the same grid and the same GLORYS T and S as the script (default inpainting), plus a
# second copy without inpainting (NaN wherever any GLORYS point used by the interpolation is
# land or below GLORYS's sea floor). No model, no sea ice, no time-stepping.
#
# Prints a census of faces with N² < -1e-6 s⁻² and saves, in ./inversions_t0/:
#   map.png         minimum N² per column, and the inpainted fraction of each column
#   site_<n>.png    at each site: horizontal slices, two vertical sections, profiles
#
# Usage, from the clone (Copernicus credentials sourced, in case a file is missing):
#   julia --project=. ../panarctic-inversions.jl                  # the three worst sites
#   julia --project=. ../panarctic-inversions.jl 179 337 148 352  # chosen (i, j) pairs
#   ARCH=CPU julia --project=. ../panarctic-inversions.jl         # without the GPU

using NumericalEarth, Oceananigans, Oceananigans.Units
using Oceananigans.ImmersedBoundaries: immersed_inactive_node
using Oceananigans.OrthogonalSphericalShellGrids: RotatedLatitudeLongitudeGrid
using Oceananigans.Grids: λnode, φnode, znode
using Oceananigans.Models: buoyancy_operation
using Dates, CUDA, Printf, Statistics
using CopernicusMarine
using CairoMakie

arch = get(ENV, "ARCH", "GPU") == "CPU" ? CPU() : GPU()

threshold     = -1e-6     # s⁻², same as the script's check
window        = 25        # half-width (cells) of the slices and sections
section_depth = 500       # m, deepest level shown in the sections
outdir        = "inversions_t0"
mkpath(outdir)

#####
##### Grid and bathymetry (as in 02_arctic_ocean_regional_script.jl)
#####

const λ₁, λ₂ = 151, 212
const φ₁, φ₂ = -34.8, 27
res = 0.15
Nx  = round(Int, (λ₂ - λ₁) / res)
Ny  = round(Int, (φ₂ - φ₁) / res)
Nz  = 20
depth = 4000meters
z = ExponentialDiscretization(Nz, -depth, 0; scale = depth/4, mutable = true)

underlying_grid = RotatedLatitudeLongitudeGrid(arch;
                                               size       = (Nx, Ny, Nz),
                                               longitude  = (λ₁, λ₂),
                                               latitude   = (φ₁, φ₂),
                                               north_pole = (145, 0),
                                               z,
                                               halo = (7, 7, 7))

dir_kw = haskey(ENV, "DATA_DIR") ? (; dir = ENV["DATA_DIR"]) : (;)

bathymetry = Metadatum(:bottom_height; dataset = ETOPO2022(), dir_kw...)
bottom_height = regrid_bathymetry(underlying_grid, bathymetry;
                                  minimum_depth = 15,
                                  interpolation_passes = 2,
                                  major_basins = 1)

grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom_height); active_cells_map = true)

#####
##### Initial T and S: as the run sees them, and without inpainting
#####

date    = DateTime(1993, 1, 1)
dataset = GLORYSDaily()
region  = BoundingBox(longitude = (-180, 180), latitude = (45, 90))

Tmd = Metadatum(:temperature; date, dataset, region, dir_kw...)
Smd = Metadatum(:salinity;    date, dataset, region, dir_kw...)

T, S       = CenterField(grid), CenterField(grid)
Traw, Sraw = CenterField(grid), CenterField(grid)

@info "Regridding GLORYS T and S (default inpainting, cached from the run)..."
set!(T, Tmd)
set!(S, Smd)
@info "Regridding GLORYS T and S without inpainting..."
set!(Traw, Tmd; inpainting = nothing)
set!(Sraw, Smd; inpainting = nothing)

# Same equation of state as ocean_simulation's default
equation_of_state = NumericalEarth.Oceans.TEOS10EquationOfState(reference_density = 1020)
buoyancy = SeawaterBuoyancy(; equation_of_state)

b    = compute!(Field(buoyancy_operation(buoyancy, grid, (; T, S))))
braw = compute!(Field(buoyancy_operation(buoyancy, grid, (T = Traw, S = Sraw))))
N²   = Array(interior(compute!(Field(∂z(b)))))    # face k lies between cells k-1 and k

T, S, Traw, Sraw, b, braw = (Array(interior(f)) for f in (T, S, Traw, Sraw, b, braw))

#####
##### Geometry on the CPU
#####

cpu_grid = Oceananigans.Architectures.on_architecture(CPU(), grid)
ccc = (Center(), Center(), Center())

wet = [!immersed_inactive_node(i, j, k, cpu_grid, ccc...) for i in 1:Nx, j in 1:Ny, k in 1:Nz]
inpainted = wet .& (isnan.(Traw) .| isnan.(Sraw))     # value at this cell came (partly) from inpainting

λ  = [λnode(i, j, 1, cpu_grid, ccc...) for i in 1:Nx, j in 1:Ny]
φ  = [φnode(i, j, 1, cpu_grid, ccc...) for i in 1:Nx, j in 1:Ny]
zc = [znode(1, 1, k, cpu_grid.underlying_grid, Center(), Center(), Center()) for k in 1:Nz]
zf = [znode(1, 1, k, cpu_grid.underlying_grid, Center(), Center(), Face())   for k in 1:Nz+1]

hb = Array(interior(bottom_height, :, :, 1))
hb[hb .≥ 0] .= NaN
nwet = dropdims(sum(wet, dims = 3), dims = 3)

#####
##### Census
#####

wet_faces = [(i, j, k) for i in 1:Nx, j in 1:Ny, k in 2:Nz if wet[i, j, k-1] && wet[i, j, k]]
faces     = sort([f for f in wet_faces if N²[f...] < threshold], by = f -> N²[f...])
from_inpainting(f) = inpainted[f...] || inpainted[f[1], f[2], f[3]-1]

println()
@printf("Wet faces checked: %d;  N² < %.0e s⁻²: %d (%.2f %%)\n",
        length(wet_faces), threshold, length(faces), 100 * length(faces) / length(wet_faces))

println("\nBy magnitude:")
edges = (-1e-6, -1e-5, -1e-4, -1e-3, -Inf)
for n in 1:4
    sel = filter(f -> edges[n+1] ≤ N²[f...] < edges[n], faces)
    nin = count(from_inpainting, sel)
    @printf("  %9.0e > N² ≥ %-6.0e  %7d faces   (%5.1f %% involve inpainted cells)\n",
            edges[n], edges[n+1], length(sel), 100 * nin / max(1, length(sel)))
end

println("\nBy level (face k between cells k-1 and k; reference depth of the face):")
for k in Nz:-1:2
    sel = filter(f -> f[3] == k, faces)
    isempty(sel) && continue
    nin = count(from_inpainting, sel)
    @printf("  k = %2d   z = %7.1f m   %7d faces   worst %.1e   %5.1f %% inpainted\n",
            k, zf[k], length(sel), N²[first(sel)...], 100 * nin / length(sel))
end

println("\nBy number of wet cells in the column:")
for (lo, hi) in ((2, 2), (3, 5), (6, 10), (11, Nz))
    sel = filter(f -> lo ≤ nwet[f[1], f[2]] ≤ hi, faces)
    @printf("  %2d–%-2d wet cells   %7d faces   (%5.1f %% inpainted)\n",
            lo, hi, length(sel), 100 * count(from_inpainting, sel) / max(1, length(sel)))
end

println("\nBy provenance:")
nin = count(from_inpainting, faces)
@printf("  one or both cells inpainted   %7d faces\n", nin)
@printf("  both cells straight GLORYS    %7d faces\n", length(faces) - nin)
@printf("\nmin S over wet cells: %.2f psu   min T: %.2f °C\n", minimum(S[wet]), minimum(T[wet]))

println("\nWhat makes the lower cell lighter:")
fresher(f) = S[f[1], f[2], f[3]-1] < S[f...]
warmer(f)  = T[f[1], f[2], f[3]-1] > T[f...]
@printf("  fresher below only   %7d\n", count(f ->  fresher(f) && !warmer(f), faces))
@printf("  warmer below only    %7d\n", count(f -> !fresher(f) &&  warmer(f), faces))
@printf("  both                 %7d\n", count(f ->  fresher(f) &&  warmer(f), faces))

println("\nClusters (20 × 20-cell boxes), the ten with most unstable faces:")
boxes = Dict{Tuple{Int, Int}, Vector{NTuple{3, Int}}}()
for f in faces
    push!(get!(boxes, (cld(f[1], 20), cld(f[2], 20)), NTuple{3, Int}[]), f)
end
for (key, sel) in first(sort(collect(boxes), by = p -> -length(p[2])), 10)
    i, j, k = first(sel)   # faces are sorted, so this is the box's worst
    @printf("  i %3d–%3d, j %3d–%3d: %6d faces, worst %.1e at (%d, %d, %d) %6.1f°E %5.1f°N, bottom %6.0f m, %3.0f %% inpainted\n",
            20key[1]-19, 20key[1], 20key[2]-19, 20key[2], length(sel), N²[i, j, k], i, j, k,
            λ[i, j], φ[i, j], hb[i, j], 100 * count(from_inpainting, sel) / length(sel))
end

println("\nThe 20 most unstable faces (below = cell k-1, above = cell k; * = inpainted):")
println("    i    j   k    lon     lat    z face  bottom      N²      T below/above      S below/above")
for f in first(faces, 20)
    i, j, k = f
    @printf("  %3d  %3d  %2d  %6.1f  %5.1f  %7.1f  %6.0f  %9.2e  %6.2f%s %6.2f%s  %6.2f%s %6.2f%s\n",
            i, j, k, λ[i, j], φ[i, j], zf[k], hb[i, j], N²[f...],
            T[i, j, k-1], inpainted[i, j, k-1] ? "*" : " ", T[f...], inpainted[f...] ? "*" : " ",
            S[i, j, k-1], inpainted[i, j, k-1] ? "*" : " ", S[f...], inpainted[f...] ? "*" : " ")
end

#####
##### Sites: (i, j) pairs from the command line, else the three worst distinct places
#####

sites = if isempty(ARGS)
    chosen = Tuple{Int, Int}[]
    for (i, j, k) in faces
        all(s -> max(abs(i - s[1]), abs(j - s[2])) > 15, chosen) && push!(chosen, (i, j))
        length(chosen) == 3 && break
    end
    chosen
else
    v = parse.(Int, ARGS)
    [(v[n], v[n+1]) for n in 1:2:length(v)-1]
end

#####
##### Map
#####

land = fill(NaN, Nx, Ny); land[isnan.(hb)] .= 1
worst_column = fill(NaN, Nx, Ny)
for (i, j, k) in faces
    worst_column[i, j] = min(isnan(worst_column[i, j]) ? 0.0 : worst_column[i, j], N²[i, j, k])
end
inpainted_fraction = [nwet[i, j] > 0 ? count(view(inpainted, i, j, :)) / nwet[i, j] : NaN for i in 1:Nx, j in 1:Ny]

fig = Figure(size = (1300, 620))
ax1 = Axis(fig[1, 1], title = "log₁₀(−N²) of the most unstable face in each column, t = 0",
           xlabel = "grid i", ylabel = "grid j", aspect = DataAspect())
heatmap!(ax1, land, colormap = [:gray75, :gray75], colorrange = (0, 2))
hm1 = heatmap!(ax1, log10.(-worst_column), colormap = :inferno, colorrange = (-6, -2))
Colorbar(fig[1, 2], hm1, label = "log₁₀(−N² / s⁻²)")
ax2 = Axis(fig[1, 3], title = "Fraction of the column's wet cells filled by inpainting",
           xlabel = "grid i", aspect = DataAspect())
heatmap!(ax2, land, colormap = [:gray75, :gray75], colorrange = (0, 2))
hm2 = heatmap!(ax2, inpainted_fraction, colormap = :viridis, colorrange = (0, 1))
Colorbar(fig[1, 4], hm2, label = "fraction")
for ax in (ax1, ax2), (n, (i, j)) in enumerate(sites)
    scatter!(ax, [i], [j], marker = :circle, color = :transparent, strokecolor = :cyan, strokewidth = 2, markersize = 18)
    text!(ax, i + 6, j + 6, text = string(n), color = :cyan, fontsize = 18)
end
save(joinpath(outdir, "map.png"), fig)

#####
##### Site figures
#####

masked(A) = (B = copy(A); B[.!wet] .= NaN; B)
Sm, Tm = masked(S), masked(T)

function site_figure(n, i₀, j₀)
    is = max(1, i₀ - window):min(Nx, i₀ + window)
    js = max(1, j₀ - window):min(Ny, j₀ + window)
    unstable = Set(faces)

    @printf("\nSite %d: (i, j) = (%d, %d), %.2f°E %.2f°N, bottom %.0f m, %d wet cells\n",
            n, i₀, j₀, λ[i₀, j₀], φ[i₀, j₀], hb[i₀, j₀], nwet[i₀, j₀])
    println("   k   z centre    T       S      T raw   S raw    N² at face k (between k-1 and k)")
    for k in Nz:-1:1
        wet[i₀, j₀, k] || continue
        n2 = k ≥ 2 && wet[i₀, j₀, k-1] ? @sprintf("%10.2e", N²[i₀, j₀, k]) : ""
        @printf("  %2d  %8.1f  %6.2f  %6.2f  %6.2f  %6.2f  %s%s\n", k, zc[k],
                T[i₀, j₀, k], S[i₀, j₀, k], Traw[i₀, j₀, k], Sraw[i₀, j₀, k], n2,
                inpainted[i₀, j₀, k] ? "   (inpainted)" : "")
    end

    window_values(A) = filter(isfinite, vec(A[is, js, max(1, Nz-8):Nz]))
    Srange = Tuple(quantile(window_values(Sm), [0.02, 0.98]))
    Trange = Tuple(quantile(window_values(Tm), [0.02, 0.98]))

    fig = Figure(size = (1700, 1300))
    Label(fig[0, 1:4], @sprintf("Site %d: (i, j) = (%d, %d), %.2f°E %.2f°N, bottom %.0f m — initial state, 1 Jan 1993",
                                n, i₀, j₀, λ[i₀, j₀], φ[i₀, j₀], hb[i₀, j₀]), fontsize = 20)

    # Row 1: horizontal slices around the site
    function slice!(position, A, title, colormap, colorrange)
        ax = Axis(fig[1, position][1, 1], title = title, xlabel = "i", ylabel = "j", aspect = DataAspect())
        heatmap!(ax, is, js, land[is, js], colormap = [:gray75, :gray75], colorrange = (0, 2))
        hm = heatmap!(ax, is, js, A, colormap = colormap, colorrange = colorrange)
        Colorbar(fig[1, position][2, 1], hm, vertical = false)
        hlines!(ax, [j₀], color = :black, linestyle = :dash)
        vlines!(ax, [i₀], color = :black, linestyle = :dash)
        scatter!(ax, [i₀], [j₀], marker = :xcross, color = :black, markersize = 14)
        return ax
    end

    ax = slice!(1, Sm[is, js, Nz], @sprintf("S, top cell (k = %d) [psu]; • inpainted", Nz), :haline, Srange)
    pts = [Point2f(i, j) for i in is, j in js if inpainted[i, j, Nz]]
    isempty(pts) || scatter!(ax, pts, color = :white, markersize = 3)

    ax = slice!(2, Sm[is, js, Nz-1], @sprintf("S, k = %d [psu]; • inpainted", Nz-1), :haline, Srange)
    pts = [Point2f(i, j) for i in is, j in js if inpainted[i, j, Nz-1]]
    isempty(pts) || scatter!(ax, pts, color = :white, markersize = 3)

    top_N² = [wet[i, j, Nz] && wet[i, j, Nz-1] ? N²[i, j, Nz] : NaN for i in is, j in js]
    N²max = max(1e-5, maximum(abs, filter(isfinite, top_N²); init = 0.0))
    slice!(3, top_N², @sprintf("N² between k = %d and %d [s⁻²]", Nz-1, Nz), :balance, (-N²max, N²max))

    slice!(4, hb[is, js], "bottom height [m]", :deep, (-min(4000, -minimum(filter(isfinite, hb[is, js]))), 0))

    # Rows 2–3: vertical sections along j = j₀ (x = i) and i = i₀ (x = j)
    function section!(row, col, A, title, colormap, colorrange, along_i)
        xs    = along_i ? is : js
        data  = along_i ? A[is, j₀, :] : A[i₀, js, :]
        floor = along_i ? hb[is, j₀] : hb[i₀, js]
        ax = Axis(fig[row, col][1, 1], title = title, xlabel = along_i ? "i" : "j", ylabel = "z [m]")
        hm = heatmap!(ax, (first(xs) - 0.5):(last(xs) + 0.5), zf, data, colormap = colormap, colorrange = colorrange,
                      nan_color = :gray75)
        Colorbar(fig[row, col][2, 1], hm, vertical = false)
        stairs!(ax, xs, floor, step = :center, color = :black)
        inp = along_i ? [Point2f(x, zc[k]) for x in xs for k in 1:Nz if inpainted[x, j₀, k]] :
                        [Point2f(x, zc[k]) for x in xs for k in 1:Nz if inpainted[i₀, x, k]]
        isempty(inp) || scatter!(ax, inp, color = :white, markersize = 4)
        bad = along_i ? [Point2f(x, zf[k]) for x in xs for k in 2:Nz if (x, j₀, k) in unstable] :
                        [Point2f(x, zf[k]) for x in xs for k in 2:Nz if (i₀, x, k) in unstable]
        isempty(bad) || scatter!(ax, bad, marker = :hline, color = :red, markersize = 14)
        vlines!(ax, [along_i ? i₀ : j₀], color = :black, linestyle = :dash)
        deepest = minimum(filter(isfinite, floor); init = 0.0)
        ylims!(ax, max(-section_depth, deepest - 10), 0)
        return ax
    end

    section!(2, 1, Sm, "S along j = $j₀ [psu]; • inpainted, red — N² < $threshold", :haline, Srange, true)
    section!(2, 2, Sm, "S along i = $i₀ [psu]", :haline, Srange, false)
    section!(3, 1, Tm, "T along j = $j₀ [°C]", :thermal, Trange, true)
    section!(3, 2, Tm, "T along i = $i₀ [°C]", :thermal, Trange, false)

    # Profiles: lines are what the run starts from; open circles are GLORYS without inpainting
    ks = [k for k in 1:Nz if wet[i₀, j₀, k]]
    function profile!(row, col, A, Araw, xlabel)
        ax = Axis(fig[row, col], xlabel = xlabel, ylabel = "z [m]")
        scatterlines!(ax, A[i₀, j₀, ks], zc[ks], color = :black, markersize = 6)
        rk = [k for k in ks if isfinite(Araw[i₀, j₀, k])]
        isempty(rk) || scatter!(ax, Araw[i₀, j₀, rk], zc[rk], color = :transparent,
                                strokecolor = :orange, strokewidth = 2, markersize = 12)
        ylims!(ax, max(-section_depth, zf[first(ks)] - 10), 0)
        return ax
    end

    profile!(2, 3, S, Sraw, "S [psu]")
    profile!(2, 4, T, Traw, "T [°C]")
    profile!(3, 3, b, braw, "b [m s⁻²]")

    ax = Axis(fig[3, 4], xlabel = "N² [s⁻²]", ylabel = "z [m]")
    kf = [k for k in 2:Nz if wet[i₀, j₀, k] && wet[i₀, j₀, k-1]]
    isempty(kf) || scatterlines!(ax, N²[i₀, j₀, kf], zf[kf], color = :black, markersize = 6)
    vlines!(ax, [0], color = :red)
    isempty(kf) || ylims!(ax, max(-section_depth, zf[first(ks)] - 10), 0)

    path = joinpath(outdir, "site_$n.png")
    save(path, fig)
    return path
end

paths = [site_figure(n, i, j) for (n, (i, j)) in enumerate(sites)]

println("\nSaved: ", join([joinpath(outdir, "map.png"); paths], ", "))
