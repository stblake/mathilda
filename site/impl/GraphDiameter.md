---
source: src/graph/gmet_distance.c
---
**Algorithm.** `builtin_graph_diameter` is the `EX_DIAMETER` case of the shared `extremal` routine. It validates the graph and reads any `EdgeWeight` list, then runs an `O(n + m)` strong-connectivity test over the out-arcs. An unweighted graph that is not strongly connected (connected, when undirected) answers `Infinity` immediately, with no all-pairs work. Otherwise it reads the cached per-source summary (`gmet_dist_summary`) and returns the largest eccentricity. Unweighted distances come from a bit-parallel multi-source BFS, 256 sources per adjacency sweep; weighted graphs use a binary-heap Dijkstra per source and a machine-real result. The empty graph gives `0`.

**Data structures.** The graph is the Expr tree `Graph[List[v1, ...], List[edge1, ...]]`; `gmet_csr_build` turns its memoized edge-index view (`graph_edge_indices`) into a CSR arc list, with an undirected edge stored as two arcs. The summary holds per-source eccentricity, reach count and distance sum arrays, so `GraphRadius`, `GraphCenter` and `MeanGraphDistance` reuse the same traversal. Results are cached by expression.

**Complexity / limits.** Unweighted cost is `O(n (n + m) / 256)` word operations, with source batches spread over the thread team; weighted cost is `O(n (m + n) log n)`. A disconnected unweighted graph is answered in `O(n + m)`. Weighted graphs return `Infinity` whenever some pair is unreachable. Symbolic, complex or negative weights leave the call unevaluated.
