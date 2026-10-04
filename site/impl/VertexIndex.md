---
source: src/graph/gops_edit.c
---
**Algorithm.** `builtin_vertex_index` answers `VertexIndex[g, v]` with the 1-based position of `v` in `VertexList[g]`, and `VertexIndex[g, {v1, v2, ...}]` with the list of positions. It calls `graph_vertex_position`, which probes the memoized vertex hash and returns `-1` for an absent vertex. A single vertex is tried first, so a vertex that is itself a list is still found as a vertex. If any item of a list argument is absent, the whole call is left unevaluated rather than returning a partial list.

**Data structures.** The graph is a `Graph[List, List]` expression tree; the vertex index is a hash from vertex expression to position, built once per graph by `graph_util.c` and reused. The result is plain integers.

**Complexity / limits.** `O(1)` expected per vertex once the index exists (`O(V)` on the first query of a graph). Returns unevaluated for a non-graph, a missing vertex or a wrong argument count.
