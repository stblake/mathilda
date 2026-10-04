---
source: src/graph/degree.c
---
**Algorithm.** `builtin_vertex_in_degree` is `degree_dispatch` in the `DEG_IN`
mode. A `DirectedEdge[a, b]` (i.e. `a -> b`) adds one to the in-degree of `b`
only; an `UndirectedEdge[a, b]` is incident to both ends and contributes one to
the in-degree of each (so in-degree equals out-degree equals total degree for a
purely undirected graph). `VertexInDegree[g, v]` returns the single integer
in-degree of `v` by the per-vertex scan `degree_of`, declining (`NULL`) when `v`
is not a vertex; `VertexInDegree[g]` returns the list of in-degrees in canonical
vertex order. A `Hypergraph` has only a total degree, so the in/out modes return
`NULL` for one.

**Data structures.** The all-vertices form makes a single `O(E)` pass that
*pushes* each edge's contribution onto its endpoints, indexing vertices through a
`GraphVIdx` open-addressing hash — this replaced an earlier pull-style
`O(V·E)` scan (~5 s on a 20000/40000 graph). Counts accumulate in an `int64`
array and are read back through the same index so a repeated vertex reports the
degree of its first occurrence; results are boxed into a `List` of integers. The
single-vertex form uses no auxiliary structure beyond the edge scan.

**Complexity / limits.** `O(V + E)` for `VertexInDegree[g]`, `O(E)` for the
single-vertex query. Self-loops are assumed absent (the constructor forbids
them), so no endpoint is double-counted within one edge. Accepts one or two
arguments; any other arity returns `NULL`.
