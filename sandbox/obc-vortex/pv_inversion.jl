# PV inversion and balanced/unbalanced split for linear rotating shallow water on a C-grid.
#
# Plain CPU arrays, no Oceananigans types: η at cell centres (Nx, Ny), u on x-faces (Nx+1, Ny), v on y-faces
# (Nx, Ny+1), H at centres (Nx, Ny). Uniform spacing Δx, Δy. Walls on the west, east and north; the south side is
# either a wall or open with a prescribed normal velocity v_ext (Nx values on the south faces).
#
# Balanced part (notes, Sections 2–3): the geostrophic state with the same linear PV q = ζ − f η / H,
#     (∇² − f²/(g H)) ψ_B = q,      η_B = f ψ_B / g,   u_B = −∂y ψ_B,   v_B = ∂x ψ_B,
# with ψ_B = ψ_w on the walls and ψ_B = ψ_w + ∫ v_ext dx on an open south side. The constant ψ_w is fixed by
# requiring the balanced part to carry the volume of the basin, ∫ η_B dA = ∫ η dA.

using SparseArrays, LinearAlgebra

struct PVInverter{F}
    Nx :: Int
    Ny :: Int
    Δx :: Float64
    Δy :: Float64
    f  :: Float64
    g  :: Float64
    H  :: Matrix{Float64}
    open_south :: Bool
    factor :: F                 # factorisation of the interior operator
    φ₁ :: Matrix{Float64}       # solution with zero PV and ψ = 1 on the whole boundary
end

# Interior operator with Dirichlet values on the cell faces of the boundary (ghost = 2b − ψ₁).
function inversion_operator(Nx, Ny, Δx, Δy, f, g, H)
    idx(i, j) = i + (j - 1) * Nx
    I = Int[]; J = Int[]; V = Float64[]
    for j in 1:Ny, i in 1:Nx
        n = idx(i, j)
        diag = -f^2 / (g * H[i, j])
        for (di, dj, Δ²) in ((-1, 0, Δx^2), (1, 0, Δx^2), (0, -1, Δy^2), (0, 1, Δy^2))
            i′, j′ = i + di, j + dj
            if 1 ≤ i′ ≤ Nx && 1 ≤ j′ ≤ Ny
                push!(I, n); push!(J, idx(i′, j′)); push!(V, 1 / Δ²)
                diag -= 1 / Δ²
            else
                diag -= 2 / Δ²          # ghost value 2b − ψ: −1/Δ² from the stencil, −1/Δ² from the ghost
            end
        end
        push!(I, n); push!(J, n); push!(V, diag)
    end
    return sparse(I, J, V, Nx * Ny, Nx * Ny)
end

# Right-hand side: q plus the boundary values (b_W, b_E, b_N scalars, b_S a vector over i).
function inversion_rhs(inv, q, bW, bE, bN, bS)
    Nx, Ny, Δx, Δy = inv.Nx, inv.Ny, inv.Δx, inv.Δy
    r = copy(q)
    r[1, :]  .-= 2bW / Δx^2
    r[Nx, :] .-= 2bE / Δx^2
    r[:, Ny] .-= 2bN / Δy^2
    r[:, 1]  .-= 2 .* bS ./ Δy^2
    return vec(r)
end

function PVInverter(Nx, Ny, Δx, Δy, f, g, H; open_south = true)
    A = inversion_operator(Nx, Ny, Δx, Δy, f, g, H)
    F = lu(A)
    tmp = PVInverter(Nx, Ny, Δx, Δy, f, g, H, open_south, F, zeros(Nx, Ny))
    φ₁ = reshape(F \ inversion_rhs(tmp, zeros(Nx, Ny), 1.0, 1.0, 1.0, ones(Nx)), Nx, Ny)
    return PVInverter(Nx, Ny, Δx, Δy, f, g, H, open_south, F, φ₁)
end

# Relative vorticity at corners (Nx+1, Ny+1), zero on the boundary corners (free slip), averaged to centres.
function vorticity_centres(u, v, Δx, Δy)
    Nx, Ny = size(u, 1) - 1, size(u, 2)
    ζ = zeros(Nx + 1, Ny + 1)
    for j in 2:Ny, i in 2:Nx
        ζ[i, j] = (v[i, j] - v[i-1, j]) / Δx - (u[i, j] - u[i, j-1]) / Δy
    end
    return [(ζ[i, j] + ζ[i+1, j] + ζ[i, j+1] + ζ[i+1, j+1]) / 4 for i in 1:Nx, j in 1:Ny]
end

potential_vorticity(inv, u, v, η) = vorticity_centres(u, v, inv.Δx, inv.Δy) .- inv.f .* η ./ inv.H

# ψ on the open south faces, relative to ψ_w: ∫₀ˣ v_ext dx at the face centres x_i.
south_streamfunction(inv, v_ext) = cumsum(v_ext .* inv.Δx) .- v_ext .* inv.Δx ./ 2

"""
    balanced_state(inv, u, v, η; v_ext = nothing)

Balanced streamfunction ψ_B (centres), the wall constant ψ_w, and the balanced fields (η_B, u_B, v_B) on the model's
C-grid. `v_ext` (south faces) is required when the south side is open.
"""
function balanced_state(inv, u, v, η; v_ext = nothing)
    Nx, Ny, Δx, Δy, f, g = inv.Nx, inv.Ny, inv.Δx, inv.Δy, inv.f, inv.g
    q  = potential_vorticity(inv, u, v, η)
    sS = inv.open_south ? south_streamfunction(inv, v_ext) : zeros(Nx)
    ψp = reshape(inv.factor \ inversion_rhs(inv, q, 0.0, 0.0, 0.0, sS), Nx, Ny)
    ψw = ((g / f) * sum(η) - sum(ψp)) / sum(inv.φ₁)          # uniform cells: area weights cancel
    ψ  = ψp .+ ψw .* inv.φ₁

    # ψ at corners: interior corners average the four centres; boundary corners take the boundary value
    sS_face = ψw .+ cumsum(vcat(0.0, (inv.open_south ? v_ext : zeros(Nx)) .* Δx))   # at x-faces 0..Nx
    ψc = zeros(Nx + 1, Ny + 1)
    for j in 2:Ny, i in 2:Nx
        ψc[i, j] = (ψ[i-1, j-1] + ψ[i, j-1] + ψ[i-1, j] + ψ[i, j]) / 4
    end
    ψc[1, :] .= ψw; ψc[Nx+1, :] .= ψw; ψc[:, Ny+1] .= ψw
    ψc[:, 1] .= sS_face
    ψc[Nx+1, 1] = ψw + (inv.open_south ? sum(v_ext) * Δx : 0.0)

    uB = [-(ψc[i, j+1] - ψc[i, j]) / Δy for i in 1:Nx+1, j in 1:Ny]
    vB = [ (ψc[i+1, j] - ψc[i, j]) / Δx for i in 1:Nx, j in 1:Ny+1]
    ηB = (f / g) .* ψ
    return (; ψ, ψw, ηB, uB, vB, sS = ψw .+ sS)
end

# Energy ½∫(H|u|² + g η²) dA with face velocities weighted by the depth at the adjacent centres.
function energy(inv, u, v, η)
    H, Δx, Δy, g = inv.H, inv.Δx, inv.Δy, inv.g
    Nx, Ny = inv.Nx, inv.Ny
    Hu(i, j) = i == 1 ? H[1, j] : i == Nx + 1 ? H[Nx, j] : (H[i-1, j] + H[i, j]) / 2
    Hv(i, j) = j == 1 ? H[i, 1] : j == Ny + 1 ? H[i, Ny] : (H[i, j-1] + H[i, j]) / 2
    Ek = sum(Hu(i, j) * u[i, j]^2 for i in 1:Nx+1, j in 1:Ny) + sum(Hv(i, j) * v[i, j]^2 for i in 1:Nx, j in 1:Ny+1)
    Ep = g * sum(abs2, η)
    return 0.5 * (Ek + Ep) * Δx * Δy
end

"""
    split_energy(inv, u, v, η; v_ext = nothing)

Total, balanced and unbalanced energy (J/ρ₀) of a state, and the balanced decomposition.
"""
function split_energy(inv, u, v, η; v_ext = nothing)
    B = balanced_state(inv, u, v, η; v_ext)
    E  = energy(inv, u, v, η)
    EB = energy(inv, B.uB, B.vB, B.ηB)
    E′ = energy(inv, u .- B.uB, v .- B.vB, η .- B.ηB)
    return (; E, EB, E′, B)
end
