---
source: src/graph/gmet_distance.c
---
**Algorithm.** `builtin_graph_periphery` gives the vertices of maximum
eccentricity — the "rim" of the graph — where a vertex's eccentricity is the
greatest distance from it to any other vertex. It shares the `extremal`
dispatcher with `GraphDiameter`/`GraphRadius`/`GraphCenter`/`MeanGraphDistance`:
an all-pairs distance summary (`gmet_dist_summary`) yields each vertex's
eccentricity, and the periphery is the set achieving the maximum. For a weighted
graph a vertex that cannot reach every other has eccentricity `Infinity` (as in
the Wolfram Language), and ties are compared with a small relative tolerance so
two equal weight-sums taken in different orders match. If `g` is not (strongly)
connected the periphery is `{}`.

**Data structures.** A `GmetCSR` adjacency plus the cached `GmetDistSummary`
(per-vertex eccentricity `ecc[]`, reach counts, distance sums). The matching
vertices are collected through `gmet_vertex_subset`, and the answer is cached on
the graph node (`gmet_cache_put`).

**Complexity / limits.** Dominated by the all-pairs distance computation (BFS per
source unweighted, Dijkstra per source weighted). Integer distances for an
unweighted graph. `{}` unless the graph is connected.
