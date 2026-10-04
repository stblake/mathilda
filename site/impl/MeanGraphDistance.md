---
source: src/graph/gmet_distance.c
---
**Algorithm.** `builtin_mean_graph_distance` is the mean distance over all ordered pairs of
distinct vertices: `tot / (n(n-1))`, where `tot` is the sum of all pairwise distances taken
from the cached per-source distance summary (`gmet_dist_summary`: multi-source BFS when
unweighted, Dijkstra when weighted). An `O(V+E)` strong-connectivity test short-circuits: a
graph that is not (strongly) connected has mean distance `Infinity`, because some pair is
unreachable.

**Data structures.** The memoised summary's per-vertex distance sums `sum[]`, shared with
`GraphCenter` and `GraphRadius`. Unweighted the result is assembled exactly as an
`Integer`/`Rational` via `make_rational((int64)tot, n(n-1))`; weighted it is a machine `Real`.

**Complexity / limits.** `O(V(V+E))` unweighted, `O(V(E + V log V))` weighted. It is undefined
for a single vertex (`n = 1`) and for the empty graph — both return unevaluated — and
`Infinity` for a non-(strongly-)connected graph; a non-graph argument returns unevaluated.
