---
source: src/graph/galg_mis.c
---
**Algorithm.** `builtin_independent_edge_set_q` tests whether a list of edges is
a **matching** of `g`: a set of edges of `g`, no two sharing a vertex. The shared
`gm_edge_set_q(res, 0)` walks the list, and for each element checks via
`gm_edge_of` that it really is an edge of `g` (accepting `DirectedEdge`/`Rule` and
`UndirectedEdge`/`TwoWayRule`, verified with `graph_has_edge`) and resolves its
two endpoints to vertex indices. A `hit[]` bitmap over the vertices records used
endpoints; if a new edge touches an already-used vertex the answer is `False`. An
element that is not an edge of `g` also gives `False`. (The same routine with the
cover flag set implements `EdgeCoverQ`.)

**Data structures.** A single `char hit[]` array of length `n` (the vertex
count) and the validated-graph membership test `graph_has_edge`; no adjacency is
built. Vertices are resolved with `galg_vertex_arg`.

**Complexity / limits.** `O(k)` edge checks for a `k`-element list, each an `O(1)`
membership probe on the graph memo. A non-graph or non-list argument gives
`False`. `FindIndependentEdgeSet` constructs a *maximum* matching.
