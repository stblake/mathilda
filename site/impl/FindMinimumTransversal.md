---
references:
  - "C. Berge, *Hypergraphs: Combinatorics of Finite Sets*, North-Holland Mathematical Library 45 (Elsevier, 1989), ch. 2 (transversals / minimum hitting set)."
source: src/graph/hyp_transversal.c
---
**Algorithm.** `builtin_find_minimum_transversal` returns one minimum-cardinality
transversal (a minimum hitting set, NP-hard) by exact **branch and bound**. It
branches on the uncovered hyperedge with the fewest open vertices — include a
vertex, then exclude it for later siblings. The search is seeded with a greedy
**max-coverage upper bound** (a lazy max-heap on coverage gains, stale entries
re-pushed, so the greedy pass is `O(Σ|e| log n)`) and pruned by
`mn_lower_bound`, the larger of three bounds: a **degree bound** (fewest
high-degree open vertices to cover all uncovered hyperedges), a **greedy disjoint
packing** (pairwise-disjoint uncovered hyperedges, taken in increasing total open
degree, one vertex each), and a **fractional-packing / LP-dual** bound
(`y_F = min_{v∈F} 1/d(v)` dual-feasible, raised by one low-degree-first
dual-ascent pass, ceiling-rounded). When a packing already matches the greedy
bound the answer is proved optimal with no search at all.

**Data structures.** The distinct-vertex CSR `soff/sv` and incidence CSR
`voff/ve`; `MinState` with the covered-count array, swap-removed uncovered set,
`excl` mask, current/best solution stacks, and the lower-bound scratch
(`odeg`/`touched`/`dcount` for the degree bound, `ekey`/`eord` ordering, `mark`
stamping for the packing, `ey`/`slack` for the LP dual). The result is a sorted
`List` of vertices.

**Complexity / limits.** Exponential worst case, kept tractable by the bounds.
**Correct-or-unevaluated:** left unevaluated if any hyperedge is empty
(unhittable), past `HYP_MIN_MAX_NODES = 2·10⁷` branch-and-bound nodes, or if the
greedy upper bound exceeds 20000 (a C-stack-depth guard); it polls
`tc_check_deadline()` every 4096 nodes for `TimeConstrained`.
