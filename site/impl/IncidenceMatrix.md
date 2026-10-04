---
source: src/graph/incmat.c
---
**Algorithm.** `builtin_incidence_matrix` builds the dense `|V| x |E|` integer incidence matrix:
row `i` is vertex `i` in `VertexList` order, column `j` is edge `j` in `EdgeList` order. For
each edge it resolves the two endpoint positions and reads the edge kind. An `UndirectedEdge
{a, b}` puts `+1` in both endpoint rows; a `DirectedEdge[a, b]` is oriented, `-1` at the tail and
`+1` at the head. Self-loops never arise, because the `Graph` constructor rejects them. For a
non-graph argument it defers to the hypergraph path `hyp_incidence_matrix`.

**Data structures.** A flat `int` grid (`calloc(n*m)`) filled column by column, then converted
to a `List` of `List`s of boxed integers — a dense matrix, not a `SparseArray`. Endpoint
positions come from a linear `graph_vertex_index` scan.

**Complexity / limits.** `O(V·E)` in both time and space (dense, plus the `O(E·V)` linear
endpoint lookups). The result is a plain nested list suitable for `MatrixForm` or the linear
algebra heads.
