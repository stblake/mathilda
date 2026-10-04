---
source: src/interval.c
---
**Algorithm.** `builtin_intervalintersection` intersects a sequence of intervals (a
bare numeric scalar is treated as a degenerate point interval). It folds left: the
accumulator starts as the pair list of the first argument (`iv_collect`, which
pushes each `[lo, hi]` pair, or `[x, x]` for a scalar/infinity, and declines on a
non-interval). For each later argument it intersects every pair of the new interval
against the whole accumulator with `iv_intersect_pair`, replacing the accumulator
with the collected overlaps. An empty accumulator yields `Interval[]` (the empty
interval — the arguments are disjoint); otherwise the surviving `[lo, hi]` pairs are
handed to `iv_canonicalize_pairs`, which sorts and merges them into the normal form.

The head is `Orderless` and `Flat`, so nested intersections flatten and argument
order does not matter before the fold runs.

**Data structures.** An `IvBuild` growable array of parallel `Expr*` lo/hi
endpoints (`ivb_init`/`ivb_push`/`ivb_free`). Endpoints are kept as exact `Expr`
leaves when exact; comparisons go through `interval_endpoint_cmp`, which can report
"undecided" for symbolic endpoints.

**Complexity / limits.** `O(p·q)` endpoint comparisons per fold step for interval
multiplicities `p`, `q`; linear in the number of arguments overall. A non-interval,
non-numeric argument makes the whole call decline (`NULL`) and stay symbolic.
Disjoint inputs correctly produce the empty `Interval[]`.
