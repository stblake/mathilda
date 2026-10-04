---
source: src/graph/degree.c
---
**Algorithm.** `builtin_vertex_out_degree` gives the out-degree of each vertex
(`VertexOutDegree[g]`, in `VertexList` order) or of a single vertex
(`VertexOutDegree[g, v]`). A `DirectedEdge[a, b]` adds `1` to `out(a)`; an
`UndirectedEdge[a, b]` is incident to both ends and so adds `1` to each of their
in-, out-, and total degrees (hence `in = out = total` for a purely undirected
graph). It shares `degree_dispatch` with `VertexDegree` and `VertexInDegree`,
selecting the `DEG_OUT` accumulator. The whole-graph form makes a single pass over
the edges, pushing each edge's contribution to its endpoints, rather than one
`O(E)` scan per vertex.

**Data structures.** For the list form, a per-vertex `int64 deg[]` accumulator and
a `GraphVIdx` hash index from vertex expression to position, so endpoints resolve
in `O(1)` and a repeated vertex reports the degree of its first occurrence. The
single-vertex form scans the edge list once (`degree_of`).

**Complexity / limits.** `O(V + E)` for the list, `O(E)` for a single vertex —
the earlier per-vertex design was `O(V*E)` (~5 s for 20000 vertices). A
`VertexOutDegree[g, v]` with `v` not a vertex leaves the call unevaluated.
