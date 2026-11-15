# Checks of the thin-bottom-cell fixes before an Oceananigans PR. One command per check:
#
#   julia --project=<PanArcticShenanigans>/nested pr_checks.jl <check> [arguments]
#
#   untouched [steps]              Grids with nothing to fix must give bit-for-bit the stock answer: flat bottom
#                                  (whole and partial), and whole cells over the plateau. "-wet" variants should be
#                                  exact on all three; the default (walls counted) changes whole cells at cliff edges.
#                                  CPU, < 1 min. Run this first.
#   mirror [days] [latitude]       The 2-D channel in x and the same channel turned into y. Every scheme must give the
#                                  same numbers both ways: tests the v-face code, which compare_schemes.jl never ran.
#   channel [days] [latitude]      compare_schemes.jl's table with the "-wet" variants added (does the fix still work
#                                  when walls do not count?).
#   seamount [days] [latitude] [gpu]  3-D flat-topped seamount: edges at every angle, staircase corners. Default
#                                  latitude 72. CPU: some minutes per run; add "gpu" to run on the GPU.
#   seam [days]                    Plateau bottom cell 20 % … 100 % wet. The fix switches on below 99 % (the 1 %
#                                  threshold): results must vary smoothly through it, with no jump at the switch.
#   nz [days]                      10, 20 and 40 levels, plateau at the same depth: the partial-cell excess should shrink
#                                  with resolution with the fix as without it.
#   static [days]                  Static levels instead of z-star.
#
# Defaults: 10 days, no rotation (seamount: 72°N). Schemes: stock (WENOVectorInvariant), hybrid, ke, flux (reference)
# and their "-wet" forms; see thin_cell_cases.jl and thickness_jump_fixes.jl. Not run in Julia before it was pushed:
# the first run of each check also tests the script.

if "gpu" in ARGS
    using CUDA
end

include(joinpath(@__DIR__, "thin_cell_cases.jl"))

argday(args, n, default) = length(args) ≥ n ? parse(Float64, args[n]) : default

#####
##### untouched: bit-for-bit
#####

function _final_state(c, bottom_kind, base, nsteps, bottom)
    model = thin_cell_model(c, bottom_kind, base; bottom)
    for _ in 1:nsteps
        time_step!(model, 10minutes)
    end
    return Dict(:u => Array(interior(model.velocities.u)), :v => Array(interior(model.velocities.v)),
                :w => Array(interior(model.velocities.w)), :b => Array(interior(model.tracers.b)),
                :η => Array(interior(model.free_surface.displacement)))
end

function final_state(c, bottom_kind, name, nsteps; bottom = bottom_function(c))
    base = use_scheme!(name)
    return Base.invokelatest(_final_state, c, bottom_kind, base, nsteps, bottom)
end

function check_untouched(args)
    nsteps = length(args) ≥ 1 ? parse(Int, args[1]) : 50
    c = ThinCellCase(geometry = :x)
    flat_on_face  = x -> -basin
    flat_partial  = x -> -2000.0          # not on a level face: every column has the same partial bottom cell
    cases = [("flat bottom, whole cells",            :whole,   flat_on_face),
             ("flat bottom, partial cells",          :partial, flat_partial),
             ("plateau, whole cells (cliff edges)",  :whole,   bottom_function(c))]
    println("Bit-for-bit against the stock scheme after $nsteps steps (largest |difference| of u, v, w, b, η).")
    println("Note: 'stock' here is WENOVectorInvariant through the override code with every switch off.\n")
    @printf("%-38s  %-11s  %-12s  %s\n", "grid", "scheme", "max |diff|", "result")
    for (label, bottom_kind, bottom) in cases
        ref = final_state(c, bottom_kind, "stock", nsteps; bottom)
        for name in ("hybrid", "hybrid-wet", "ke", "ke-wet")
            s = final_state(c, bottom_kind, name, nsteps; bottom)
            d = maximum(maximum(abs, s[f] .- ref[f]) for f in keys(ref))
            expected = startswith(name, "ke") ? "(ke includes symdiv everywhere: change expected)" :
                       name == "hybrid" && bottom_kind == :whole && label[1] == 'p' ? "(walls counted: change expected)" : ""
            @printf("%-38s  %-11s  %-12.3e  %s %s\n", label, name, d, d == 0 ? "identical" : "differs", expected)
        end
    end
    println("""

Pass: "hybrid-wet" identical on all three grids; "hybrid" identical on the two flat bottoms. If so, the PR can say
that grids without thin-versus-thick bottom cells are unchanged bit for bit (with the wet-cells-only detection).""")
end

#####
##### mirror: x versus y
#####

function check_mirror(args)
    days = argday(args, 1, 10.0)
    latitude = argday(args, 2, 0.0)
    schemes = ["stock", "hybrid", "hybrid-wet", "ke"]
    cx = ThinCellCase(geometry = :x, latitude = latitude)
    cy = ThinCellCase(geometry = :y, latitude = latitude)
    println("Mirror test: the same channel along x and along y, $days days, ",
            latitude == 0 ? "no rotation" : "f at $(latitude)°N", ".\n")
    @printf("%-11s %-8s  %-12s %-12s %-10s  %-9s %-9s  %-10s %-10s  %s\n", "scheme", "bottom", "mixing x", "mixing y",
            "rel. diff", "jet x", "jet y", "energy x", "energy y", "result")
    worst = 0.0
    for name in schemes, bottom_kind in (:whole, :partial)
        rx = run_case(cx, bottom_kind, name; days)
        ry = run_case(cy, bottom_kind, name; days)
        rel = max(abs(rx.mixing - ry.mixing) / abs(rx.mixing), abs(rx.jet - ry.jet) / rx.jet,
                  abs(rx.energy - ry.energy) / max(abs(rx.energy), 1e-12))
        worst = max(worst, rel)
        @printf("%-11s %-8s  %-12.6f %-12.6f %-10.1e  %-9.5f %-9.5f  %+-10.5f %+-10.5f  %s\n", name, bottom_kind,
                rx.mixing, ry.mixing, rel, rx.jet, ry.jet, rx.energy, ry.energy, rel < 1e-3 ? "same" : "DIFFERENT")
    end
    println("""

Pass: every row "same" (relative difference below 1e-3; round-off alone gives ~1e-12 or less). A difference only for
hybrid/ke partial points at the v-face code of the fix (thickness_changes_y, flux_form_v∂yv, U_dot_∇v).""")
    @printf("Largest relative difference: %.1e\n", worst)
end

#####
##### channel, seamount, nz, static: scheme tables
#####

check_channel(args) = scheme_table(ThinCellCase(geometry = :x, latitude = argday(args, 2, 0.0)),
                                   ["stock", "hybrid", "hybrid-wet", "ke", "ke-wet", "flux"]; days = argday(args, 1, 10.0))

function check_seamount(args)
    arch = "gpu" in args ? GPU() : CPU()
    numbers = filter(a -> a != "gpu", args)
    c = ThinCellCase(geometry = :xy, latitude = argday(numbers, 2, 72.0), arch = arch)
    scheme_table(c, ["stock", "hybrid", "hybrid-wet", "ke", "flux"]; days = argday(numbers, 1, 10.0))
    grid = thin_cell_grid(with(c; arch = CPU()), :partial)
    for walls in (true, false)
        set_jump_walls!(walls)
        print(walls ? "walls counted:    " : "wet cells only:   ")
        Base.invokelatest(report_fix_faces, grid)
    end
end

function check_seam(args)
    days = argday(args, 1, 10.0)
    c = ThinCellCase(geometry = :x)
    φs = [0.2, 0.5, 0.9, 0.97, 0.985, 0.995, 1.0]
    schemes = ["stock", "hybrid", "hybrid-wet"]
    println("Seam test: plateau bottom cell from 20 % to 100 % wet, partial cells, $days days. The fix acts at the plateau")
    println("level when the cell is less than 99 % wet (1 % threshold); between 0.985 and 0.995 it switches off there.\n")
    whole = run_case(c, :whole, "stock"; days)
    @printf("whole cells, stock (reference): mixing %.4f, jet %.4f m/s\n\n", whole.mixing, whole.jet)
    @printf("%-7s %-9s", "φ", "plateau")
    for s in schemes
        @printf("  %-22s", "$s: mixing / jet")
    end
    println()
    rows = []
    for φ in φs
        cφ = with(c; plateau = plateau_for_fraction(φ, c))
        r = [run_case(cφ, :partial, s; days) for s in schemes]
        push!(rows, r)
        @printf("%-7.3f %-9.1f", wet_fraction(cφ), cφ.plateau)
        for x in r
            @printf("  %8.4f / %-11.4f", x.mixing, x.jet)
        end
        println()
    end
    println("""

Pass: each column changes smoothly with φ, and hybrid's step between φ = 0.985 and 0.995 is no larger than stock's
between the same rows (at φ ≥ 0.99 the cell is only 1 % thinner than its neighbour, so the switch should hardly matter).""")
    i₁, i₂ = findfirst(==(0.985), φs), findfirst(==(0.995), φs)
    for (n, s) in enumerate(schemes)
        @printf("  %-11s step across the switch: mixing %+.4f, jet %+.4f m/s\n", s,
                rows[i₂][n].mixing - rows[i₁][n].mixing, rows[i₂][n].jet - rows[i₁][n].jet)
    end
end

function check_nz(args)
    days = argday(args, 1, 10.0)
    schemes = ["stock", "hybrid", "hybrid-wet", "flux"]
    println("Resolution test: plateau at $(default_plateau) m with 10, 20 and 40 levels (the wet fraction changes with Nz).")
    for Nz in (10, 20, 40)
        scheme_table(ThinCellCase(geometry = :x, Nz = Nz), schemes; days)
    end
    println("""

Pass: with the fix the partial/whole ratio stays near 1 at every Nz, and partial-cell mixing does not grow with Nz.""")
end

check_static(args) = scheme_table(ThinCellCase(geometry = :x, vertical = :static, latitude = argday(args, 2, 0.0)),
                                  ["stock", "hybrid", "hybrid-wet", "flux"]; days = argday(args, 1, 10.0))

#####
##### Dispatch
#####

const checks = Dict("untouched" => check_untouched, "mirror" => check_mirror, "channel" => check_channel,
                    "seamount" => check_seamount, "seam" => check_seam, "nz" => check_nz, "static" => check_static)

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("usage: pr_checks.jl <$(join(sort(collect(keys(checks))), "|"))> [arguments]")
    haskey(checks, ARGS[1]) || error("unknown check $(ARGS[1]); one of $(join(sort(collect(keys(checks))), ", "))")
    t = @elapsed checks[ARGS[1]](ARGS[2:end])
    @printf("\n(%s took %.1f min)\n", ARGS[1], t / 60)
end
