---
source: src/graph/gops_edit.c
---
**Algorithm.** `builtin_vertex_add` takes `VertexAdd[g, v]` or `VertexAdd[g, {v1, v2, ...}]`. It first filters the items: any item already a vertex of `g` (a lookup in the memoized vertex index) is dropped, and repeats among the new items are dropped through a scratch hash, so the survivors keep their first-appearance order. If nothing is new, a copy of `g` is returned. Otherwise the new vertices are appended after the existing ones; the edge list, endpoint arrays and any `EdgeWeight` list are carried over unchanged.

**Data structures.** `Graph[List[vertices], List[edges]]` expression tree; the edit reads the memoized endpoint arrays (`gops_view`) and constructs the result with `gops_graph_new`. Vertex and edge nodes are copied by reference count, not deep-copied.

**Complexity / limits.** `O(V + E)` for the copy plus one hash probe per argument item. Added vertices are isolated, so they cannot create loops or parallel edges. A first argument that is not a valid graph, or a call with other than two arguments, is left unevaluated.
