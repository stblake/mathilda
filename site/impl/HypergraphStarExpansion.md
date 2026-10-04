---
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hypergraph_star_expansion` gives the incidence (star)
bipartite `Graph`. Its vertices are `VertexList[h]` followed by one node
`Hyperedge[j]` per hyperedge; for each vertex `v` in the distinct-vertex set of
hyperedge `j` it emits `v <-> Hyperedge[j]`. Before building, each `Hyperedge[j]`
node is looked up in the hypergraph's own vertex index; if some vertex of `h` *is*
such a node the construction would be ambiguous, so the call is left unevaluated
(`clash`). The argument may be a `Hypergraph` or a bare List of hyperedges
(`hyp_arg`).

**Data structures.** The distinct-vertex CSR `soff/sv` and the memo's vertex
`GraphVIdx` (for the clash check). Vertices and edges are gathered into `Expr**`
arrays; the result is a `Graph`. Edge count equals the total incidence
`Σ|distinct e_j|`.

**Complexity / limits.** `O(n + Σ|e|)` — one `Hyperedge[j]` node per hyperedge and
one edge per incidence. The clash guard is the only reason the head declines a
well-formed hypergraph.
