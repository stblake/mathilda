---
references:
  - "P. Hage and F. Harary, *Eccentricity and centrality in networks*, Social Networks **17** (1995) 57-63."
source: src/graph/gmet_centrality.c
---
**Algorithm.** `builtin_eccentricity_centrality` is a one-line call to
`summary_centrality(res, 1)`. That helper fetches the per-source distance summary
`gmet_dist_summary(g)` and returns, for every vertex `v`, the reciprocal
eccentricity `1/ecc[v]` (and `0` when `ecc[v] = 0`, i.e. an isolated vertex), where
`ecc[v]` is the largest shortest-path distance from `v` to any vertex reachable
along its out-arcs. The summary itself does all the work: it builds an out-arc CSR
and, for each source, records the number of reachable vertices, the sum of their
distances, and their maximum. Unweighted graphs use a bit-parallel multi-source
BFS (a block of sources advanced together with a bitset frontier, parallelised over
blocks); when an `EdgeWeight` is present it switches to a binary-heap Dijkstra per
source, using the weights as lengths.

**Data structures.** A `GmetCSR` compressed-sparse-row adjacency, and a
`GmetDistSummary` (`reach[]`, `sum[]`, `ecc[]`) cached per graph node in a small
slot table and evicted round-robin, so a second centrality query on the same graph
reuses the distances. The result is handed back through `gmet_real_vector`, which
packs into an `f64` NDArray buffer above the packing threshold and otherwise builds
a `List` of reals.

**Complexity / limits.** `O(V(V+E))` unweighted (the BFS amortised over the
word-wide source blocks) and `O(V · E log V)` weighted, both multithreaded over the
source set. The value is a machine real, never exact. A non-graph argument, or one
carrying unusable weights, leaves the call unevaluated (`NULL`). In a disconnected
graph the eccentricity is taken over the reachable set only, so a vertex's score
reflects its own component.
