---
source: src/graph/vertexlist.c
---
**Algorithm.** `builtin_vertex_list` is a thin reader: for a valid graph it returns a copy of
the graph's vertex `List` — the first argument of the canonical `Graph[List[verts],
List[edges]]` node — in canonical order, with no computation. For a non-graph argument it defers
to the hypergraph reader `hyp_vertex_list`, which handles a `Hypergraph[...]` argument or returns
unevaluated.

**Data structures.** None beyond the deep copy of the stored vertex list; the canonical graph
already holds its vertices in the order `VertexList` reports.

**Complexity / limits.** `O(V)` for the copy. Always exactly one argument.
