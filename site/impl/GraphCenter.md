---
source: src/graph/gmet_distance.c
---
**Algorithm.** `builtin_graph_center` gives the vertices whose eccentricity equals the graph
radius (the minimum eccentricity). Through the shared `extremal_compute`, it reads the cached
per-source distance summary (`gmet_dist_summary`: multi-source BFS unweighted, Dijkstra
weighted), finds the minimum eccentricity `target`, and returns every vertex whose eccentricity
matches it (within a relative tolerance `ECC_TIE_TOL = 1e-12` on the weighted path). An
`O(V+E)` strong-connectivity test short-circuits first: for an unweighted graph that is not
(strongly) connected the centre is empty.

**Data structures.** The memoised distance summary's `ecc[]` array, shared with `GraphRadius`
and `MeanGraphDistance`; `gmet_vertex_subset` lifts the flagged indices back to the vertex
expressions in `VertexList` order.

**Complexity / limits.** `O(V(V+E))` unweighted, `O(V(E + V log V))` weighted. The empty graph
and a non-(strongly-)connected unweighted graph both return `{}`; a non-graph argument returns
unevaluated.
