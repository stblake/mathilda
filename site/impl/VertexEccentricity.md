---
source: src/graph/gmet_distance.c
---
**Algorithm.** `builtin_vertex_eccentricity` gives `VertexEccentricity[g, v]`, the
largest distance from `v` to any other vertex. It runs a single-source distance
computation from `v` (`single_source` — BFS for an unweighted graph, Dijkstra
when `g` carries usable `EdgeWeight`) and takes the maximum. For an **unweighted**
graph the maximum is over the vertices `v` can reach; for a **weighted** graph
Wolfram measures over *all* vertices, so a vertex `v` cannot reach makes the
eccentricity `Infinity`. The result is an `Integer` for an unweighted graph and a
machine `Real` otherwise.

**Data structures.** The CSR distance machinery of `src/graph/gmet_distance.c`:
an `int64 di[]` array of hop distances and a `double dw[]` of weighted distances,
filled by the single-source routine from the memo's integer endpoints and the
resolved edge weights (`gmet_edge_weights`).

**Complexity / limits.** `O(V + E)` unweighted, `O(E + V log V)` weighted, per
call. The graph-wide maxima/minima of eccentricity are `GraphDiameter`,
`GraphRadius`, `GraphCenter`, and `GraphPeriphery`, which compute all
eccentricities at once and cache them. A `v` that is not a vertex leaves the call
unevaluated.
