---
source: src/graph/gmet_distance.c
---
**Algorithm.** `builtin_graph_radius` is the minimum vertex eccentricity. Through the shared
`extremal_compute`, it reads a cached per-source distance summary (`gmet_dist_summary`:
multi-source BFS when unweighted, Dijkstra when weighted), computes each vertex's eccentricity
`ecc[i]` (the largest distance from `i` to a vertex it reaches), and returns the minimum. An
`O(V+E)` strong-connectivity test short-circuits first: an unweighted graph that is not
(strongly) connected has radius `Infinity`. A weighted graph can still return a finite radius
even when some pair is unreachable, following Wolfram's rule that a weighted eccentricity is
`Infinity` only when the vertex fails to reach all `n-1` others.

**Data structures.** The memoised distance summary's `ecc[]`/`reach[]` arrays, shared with
`GraphCenter` and `MeanGraphDistance` so the three heads cost one traversal between them. The
value is returned exactly (an `Integer`) for an unweighted graph and as a machine `Real` for a
weighted one, through `gmet_distance_value`. A weighted eccentricity tie uses a relative
tolerance `ECC_TIE_TOL = 1e-12`.

**Complexity / limits.** `O(V(V+E))` unweighted, `O(V(E + V log V))` weighted. The empty graph
returns `0`; a non-(strongly-)connected unweighted graph returns `Infinity`; a non-graph
argument returns unevaluated.
