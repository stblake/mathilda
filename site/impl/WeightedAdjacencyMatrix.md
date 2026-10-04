---
source: src/graph/wtadjmat.c
---
**Algorithm.** `builtin_weighted_adjacency_matrix` builds the dense `n x n`
adjacency matrix filled with per-edge **weights** instead of a literal `1`. It is
the same construction as `AdjacencyMatrix`: a `DirectedEdge[a, b]` sets
`M[a][b] = weight(a, b)`, an `UndirectedEdge` sets both `M[a][b]` and `M[b][a]`,
and any cell with no edge is `0`. The weights come from
`graph_resolve_edge_weights`, the shared resolver that returns `g`'s `EdgeWeight`
list or `{1, ..., 1}` when there is none — so for an unweighted graph
`WeightedAdjacencyMatrix[g]` equals `AdjacencyMatrix[g]` exactly.

**Data structures.** A flat `n*n` grid of `Expr*` cells, a `GraphVIdx` mapping
each vertex to its row/column index, and the resolved weight list. The grid is
assembled into a `List` of row `List`s; empty cells become `expr_new_integer(0)`.

**Complexity / limits.** `O(n^2)` to materialize the dense matrix plus `O(E)` to
place the weights. Weights are copied verbatim, so symbolic weights appear as-is
in the matrix. Returns `NULL` (unevaluated) for a non-graph argument.
