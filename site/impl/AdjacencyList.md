---
source: src/graph/adjlist.c
---
**Algorithm.** `builtin_adjacency_list` returns the neighbours of each vertex, following
the same orientation convention as `AdjacencyMatrix`: for a `DirectedEdge[v, u]` only the
successor `u` is a neighbour of `v`, while an `UndirectedEdge` contributes both endpoints.
`AdjacencyList[g]` gives `{neighbours(v1), neighbours(v2), ...}` in canonical vertex order;
`AdjacencyList[g, v]` gives just the neighbours of `v`. For each vertex the helper
`neighbors_of` scans `g`'s edge list once, comparing endpoints by `expr_eq` and appending
each new neighbour in **first-appearance order**, de-duplicated (so a vertex reachable by
both an undirected and a redundant directed edge is listed once).

**Data structures.** The graph is the canonical `Graph[List verts, List edges]` `Expr`
tree, read directly — `neighbors_of` allocates a scratch `Expr*` buffer of size `2·|E|`,
`expr_copy`s each surviving neighbour into a fresh `List`, and the one-argument form wraps
the per-vertex lists in an outer `List`. Vertex membership for the two-argument form is
checked with `graph_vertex_index` (a linear `expr_eq` scan over `verts`).

**Complexity / limits.** `neighbors_of` is `O(E)` per vertex with an inner dedup scan, so
`AdjacencyList[g]` is `O(V·E)` in the worst case (it reads the raw edge list rather than
the shared `GraphAdj` CSR). The head returns `NULL` (unevaluated) when `g` is not a valid
graph, and `AdjacencyList[g, v]` returns `NULL` when `v` is not a vertex of `g`.
