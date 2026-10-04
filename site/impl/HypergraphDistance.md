---
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hypergraph_distance` gives the least number of hyperedges
in a chain of vertices `u = x_0, x_1, ..., x_k = v` with consecutive vertices
sharing a hyperedge — the shortest-path distance in the clique expansion
(`HypergraphCliqueExpansion`), `Infinity` when none, `0` for `u == v`. It runs a
breadth-first search over vertices that, instead of materialising the clique
graph, expands each hyperedge at most once (an `eexp` stamp): popping a vertex
walks its incidence list, and each not-yet-expanded hyperedge relaxes all its
distinct members. `HypergraphDistance[h, u]` returns distances from `u` to every
vertex in `VertexList` order.

**Data structures.** The distinct-vertex CSR `soff/sv` and incidence CSR
`voff/ve`; a BFS queue, a vertex distance array, and an `eexp` hyperedge-expanded
bitmap; `dist_list`/`dist_expr` render the result (packed `int64` when all finite,
`Infinity` symbols otherwise).

**Complexity / limits.** `O(Σ|e|)` — each hyperedge's members are scanned once.
An unknown source or target vertex leaves the call unevaluated.
