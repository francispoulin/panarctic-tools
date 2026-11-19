Thanks Simone, nice idea. I ran our thin-bottom-cell cases (from #6532/#6533) on this branch and on its base (acfcbdf), 10 days, partial vs whole cells, plateau bottom cell 20 % wet. Partial/whole mixing ratio and near-bottom jet:

| Case | main | #6538 | #6533 | flux form |
|---|---|---|---|---|
| 2-D channel, f = 0 | 2.39, 0.40 m/s | 1.24, 0.25 | 1.15, 0.25 | 1.12, 0.27 |
| 2-D channel, 72°N (x and y identical) | 1.15, 0.37 | 1.06, 0.28 | 1.04, 0.23 | 1.02, 0.23 |
| 3-D seamount, 72°N | 2.11, 0.85 | 1.13, 0.39 | 1.03, 0.27 | 1.02, 0.30 |

Whole cells are unchanged in every case and the energy now falls everywhere, so this removes most of the problem while keeping one formulation, which I agree is better. I'm happy to approve it, and I'll close #6533 in favour of this PR.

What remains is clearest on the seamount. When we diagnosed #6532 we found two inconsistencies at a thickness jump: the divergence in the vertical term (which this PR fixes) and the kinetic-energy gradient, which averages u² over faces of very different thickness. I suspect the residual is the second. Would a similar treatment of the kinetic-energy gradient at thickness jumps fit into this PR, or a follow-up?

The cases run with plain Oceananigans on any branch: [`branch_checks.jl`](https://github.com/francispoulin/panarctic-tools/blob/main/sandbox/lake_at_rest/branch_checks.jl) (CPU, a few minutes).
