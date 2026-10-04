---
source: src/graph/gops_edit.c
---
**Algorithm.** `builtin_edge_index` gives the 1-based position of an edge in
`EdgeList[g]`. A single edge argument — written `u -> v`, `u <-> v`, or an
explicit `DirectedEdge`/`UndirectedEdge` — is resolved to its `EdgeList`
position; an undirected edge matches either orientation. `EdgeIndex[g, {e1,
...}]` maps the resolver over a list and returns the list of positions. Any edge
not present (or a mixed-up orientation for a directed edge) leaves the call
unevaluated.

**Data structures.** It opens a `GopsView` over the graph (the shared borrowed
view of vertices, integer endpoint arrays `eu[]`/`ev[]`, and direction flags),
parses each query edge to endpoint positions with `gops_parse_edge`, and resolves
them through the subsystem's `resolve_edges` helper, which indexes `g`'s edges so
each lookup is a hash probe rather than a scan.

**Complexity / limits.** `O(V + E)` to build the view and edge index, then `O(1)`
per queried edge. The companion `VertexIndex` does the same for vertices. Returns
`NULL` (unevaluated) if any queried edge is absent.
