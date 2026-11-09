"""Check of section_monitor.jl's path logic on the model grid geometry (polar stereographic, as ConformalConicGrid builds it).

1. Corner lon/lat from Oceananigans' formulas; the old jet band's ends.
2. The same corner_path / path_faces algorithm; a discrete streamfunction ψ on corners gives face velocities whose
   section transport must equal ψ(end) − ψ(start) exactly (path connected, signs right), and a back-and-forth
   staircase must not change in/out parts.
"""
import numpy as np

R = 6.371e6
λ0 = -40.0
spacing = 15e3
x1, x2, y1, y2 = -3280e3, 3616e3, -4100e3, 2700e3
Nx, Ny = round((x2 - x1) / spacing), round((y2 - y1) / spacing)
Δx, Δy = (x2 - x1) / Nx, (y2 - y1) / Ny
print(f"grid {Nx} x {Ny}, Δx {Δx/1e3:.3f} km, Δy {Δy/1e3:.3f} km")

def to_geo(x, y):
    ρ = np.hypot(x, y)
    θ = np.arctan2(x, -y)
    φ = 2 * np.arctan(2 * R / ρ) - np.pi / 2       # T = (R F / ρ)^(1/n), φ = 2 atan T − π/2, F = 2, n = 1
    λ = np.degrees(θ) + λ0
    return (λ + 180) % 360 - 180, np.degrees(φ)

# corners (i, j), 1-based: x = x1 + (i - 1) Δx
I, J = np.meshgrid(np.arange(1, Nx + 2), np.arange(1, Ny + 2), indexing="ij")
cλ, cφ = to_geo(x1 + (I - 1) * Δx, y1 + (J - 1) * Δy)
def corner(i, j): return cλ[i - 1, j - 1], cφ[i - 1, j - 1]

print("old jet band ends: corner (335, 175) = %.2f°E %.2f°N, corner (353, 175) = %.2f°E %.2f°N" % (*corner(335, 175), *corner(353, 175)))
print("Atlantic edge corners: %.1f°E %.1f°N ... %.1f°E %.1f°N" % (*corner(1, 1), *corner(Nx + 1, 1)))

def gc(λ1, φ1, λ2, φ2):
    r = np.radians
    return np.degrees(2 * np.arcsin(np.sqrt(np.sin(r(φ2 - φ1) / 2) ** 2 +
                                            np.cos(r(φ1)) * np.cos(r(φ2)) * np.sin(r(λ2 - λ1) / 2) ** 2)))

def nearest(λ, φ, near=None, reach=4):
    if near is None:
        d = gc(λ, φ, cλ, cφ)
    else:
        i0, j0 = near
        ia, ib = max(1, i0 - reach), min(Nx + 1, i0 + reach)
        ja, jb = max(1, j0 - reach), min(Ny + 1, j0 + reach)
        d = np.full(cλ.shape, np.inf)
        d[ia-1:ib, ja-1:jb] = gc(λ, φ, cλ[ia-1:ib, ja-1:jb], cφ[ia-1:ib, ja-1:jb])
    k = np.unravel_index(np.argmin(d), d.shape)
    return (k[0] + 1, k[1] + 1)

grid_spacing = min(gc(*corner(i, j), *corner(i + 1, j)) for i in range(1, Nx + 1) for j in (1, Ny // 2, Ny + 1))

def corner_path(points):
    path = [nearest(*points[0])]
    for (λa, φa), (λb, φb) in zip(points[:-1], points[1:]):
        n = int(np.ceil(8 * gc(λa, φa, λb, φb) / grid_spacing))
        for s in range(1, n + 1):
            α = s / n
            c = nearest((1 - α) * λa + α * λb, (1 - α) * φa + α * φb, near=path[-1])
            while c != path[-1]:
                i, j = path[-1]
                path.append((i + int(np.sign(c[0] - i)), j) if i != c[0] else (i, j + int(np.sign(c[1] - j))))
    return path

def path_faces(path):
    signs = {}
    for a, b in zip(path[:-1], path[1:]):
        Δi, Δj = b[0] - a[0], b[1] - a[1]
        assert abs(Δi) + abs(Δj) == 1, "path not made of unit steps"
        face = ("v", min(a[0], b[0]), a[1]) if Δj == 0 else ("u", a[0], min(a[1], b[1]))
        signs[face] = signs.get(face, 0) + (Δi if Δj == 0 else -Δj)
    return [(*f, s) for f, s in signs.items() if s != 0]

# Discrete streamfunction on corners; face transports (per unit depth) consistent with "transport to the left of
# a segment a → b = ψ(b) − ψ(a)":
#   v-face (i, j) (segment (i,j)→(i+1,j), left = +j):  v = ψ(i+1, j) − ψ(i, j)
#   u-face (i, j) (segment (i,j)→(i,j+1), left = −i):  −u = ψ(i, j+1) − ψ(i, j)
rng = np.random.default_rng(1)
ψ = rng.standard_normal((Nx + 1, Ny + 1)).cumsum(0).cumsum(1)
v = ψ[1:, :] - ψ[:-1, :]          # (Nx, Ny+1), v[i-1, j-1] is v-face (i, j)
u = -(ψ[:, 1:] - ψ[:, :-1])       # (Nx+1, Ny), u[i-1, j-1] is u-face (i, j)

def transport(faces):
    q = [s * (v[i - 1, j - 1] if k == "v" else u[i - 1, j - 1]) for k, i, j, s in faces]
    return sum(q), sum(x for x in q if x > 0), sum(x for x in q if x < 0)

sections = [
    ("fram",    [(-22.0, 78.83), (13.0, 78.83)]),
    ("bso",     [(16.5, 76.6), (19.0, 74.45), (20.6, 70.1)]),
    ("svinoy",  [(-1.0, 65.0), (6.2, 62.25)]),
    ("lofoten", [(9.25, 69.71), (13.31, 67.88)]),
    ("davis",   [(-63.0, 66.6), (-52.5, 66.6)]),
]
for name, pts in sections:
    path = corner_path(pts)
    faces = path_faces(path)
    net = transport(faces)[0]
    exact = ψ[path[-1][0] - 1, path[-1][1] - 1] - ψ[path[0][0] - 1, path[0][1] - 1]
    length = sum(gc(*corner(*a), *corner(*b)) for a, b in zip(path[:-1], path[1:])) * np.pi / 180 * R / 1e3
    straight = sum(gc(*a, *b) for a, b in zip(pts[:-1], pts[1:])) * np.pi / 180 * R / 1e3
    print(f"{name:8s} corners {path[0]} → {path[-1]}, {len(faces):3d} faces, staircase {length:6.0f} km "
          f"(line {straight:5.0f} km), ψ check {net - exact:+.1e}, ends "
          f"{corner(*path[0])[0]:.2f}°E {corner(*path[0])[1]:.2f}°N / {corner(*path[-1])[0]:.2f}°E {corner(*path[-1])[1]:.2f}°N")

# old band: row j = 175 corners 335 → 353 must be exactly the v-faces i = 335..352, j = 175, sign +1
path = corner_path(sections[3][1])
faces = sorted(path_faces(path))
old = [("v", i, 175, 1) for i in range(335, 353)]
print("lofoten section = old jet band faces:", faces == old, "" if faces == old else faces[:5])

# back-and-forth: a path that steps across a face and back must not count it
f = path_faces([(10, 10), (11, 10), (11, 11), (11, 10), (12, 10)])
print("back-and-forth faces:", f)
