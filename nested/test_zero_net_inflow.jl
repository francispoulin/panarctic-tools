# The branch's test of zero_net_inflow ("zero_net_inflow removes the net inflow of the Flather exterior", in
# NumericalEarth.jl/test/test_nested_ocean_model.jl), runnable in this environment: NumericalEarth's own test
# environment pulls in Breeze, SpeedyWeather, Python and GDAL, which we do not need. Keep the two in step.
#
#   julia --project=. test_zero_net_inflow.jl          (CPU)
#   GPU_TEST=true julia --project=. test_zero_net_inflow.jl

using NumericalEarth
using Oceananigans
using Oceananigans.Architectures: on_architecture
using Oceananigans.BoundaryConditions: getbc
using Oceananigans.Operators: Δxᶜᶠᶜ, Δyᶠᶜᶜ
using Oceananigans.Units: days
using CUDA
using Test

arch = get(ENV, "GPU_TEST", "false") == "true" ? GPU() : CPU()

function nesting_test_parent(arch; u = 0.1, T = 15, S = 35, η = 0)
    grid = LatitudeLongitudeGrid(arch; size = (12, 12, 6),
                                 longitude = (-2, 2), latitude = (28, 32), z = (-1000, 0),
                                 topology = (Bounded, Bounded, Bounded))
    parent_ocean = PrescribedOcean(grid, [0.0, 10days])
    set!(parent_ocean; T, S, u, v = 0, η)
    return parent_ocean
end

nesting_test_child_grid(arch) = LatitudeLongitudeGrid(arch; size = (20, 20, 6), halo = (7, 7, 7),
                                                      longitude = (-1, 1), latitude = (29, 31), z = (-1000, 0),
                                                      topology = (Bounded, Bounded, Bounded))

@testset "zero_net_inflow removes the net inflow of the Flather exterior on $(arch)" begin
    # The parent's u = 0.1 m/s enters through the west side and nothing leaves through the north side
    parent_ocean = nesting_test_parent(arch)
    parent_bottom_height = Field{Center, Center, Nothing}(parent_ocean.grid)
    set!(parent_bottom_height, -1000)

    function net_exterior_inflow(zero_net_inflow)
        model = nested_ocean_model(parent_ocean, nesting_test_child_grid(arch); parent_bottom_height,
                                   open_sides = (:west, :north), zero_net_inflow)
        child = model.child
        grid = on_architecture(CPU(), child.grid)
        Nx, Ny, _ = size(grid)
        barotropic_velocities = child.free_surface.barotropic_velocities
        U = on_architecture(CPU(), barotropic_velocities.U.boundary_conditions.west.condition)
        V = on_architecture(CPU(), barotropic_velocities.V.boundary_conditions.north.condition)
        west = sum(first(getbc(U, j, 1, grid, child.clock)) * Δyᶠᶜᶜ(1, j, 1, grid) for j in 1:Ny)
        north = sum(first(getbc(V, i, 1, grid, child.clock)) * Δxᶜᶠᶜ(i, Ny + 1, 1, grid) for i in 1:Nx)
        return west - north
    end

    inflow = net_exterior_inflow(false)
    @info "Net exterior inflow: $(inflow / 1e6) Sv without, $(net_exterior_inflow(true) / 1e6) Sv with zero_net_inflow"
    @test inflow > 0
    @test abs(net_exterior_inflow(true)) < 1e-10 * inflow
end
