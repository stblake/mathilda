---
source: src/graph/degree.c
---
**Algorithm.** `degree_dispatch` serves `VertexDegree[g]` and `VertexDegree[g, v]`, and shares its code with `VertexInDegree` and `VertexOutDegree`. A `DirectedEdge[a, b]` adds 1 to the out-degree of `a` and the in-degree of `b`, and total degree is in plus out. An `UndirectedEdge[a, b]` adds 1 to each endpoint's in, out and total degree. There are no self-loops, so no edge is counted twice at one vertex. A hypergraph argument is handed to `hyp_vertex_degree`.

**Data structures.** The list form makes a single pass over the edges. It resolves endpoints through a `GraphVIdx` hash index into an `int64_t` counter array, then emits a `List` of integers in canonical `VertexList` order. This replaced one `O(E)` scan per vertex, which took about 5 s at 20000 vertices and 40000 edges. The single-vertex form does one `expr_eq` scan over the edge list.

**Complexity / limits.** `O(V + E)` for the full list and `O(E)` for one vertex. The result is not cached. If `v` is not a vertex of `g`, or the argument is not a graph or hypergraph, the call is left unevaluated.
