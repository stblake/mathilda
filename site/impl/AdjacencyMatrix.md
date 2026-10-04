---
source: src/graph/adjmat.c
---
**Algorithm.** `builtin_adjacency_matrix` takes exactly one argument and declines (`NULL`) unless `graph_is_valid` accepts it. It allocates an `n x n` zeroed `int` grid, resolves each edge endpoint to a vertex position, and sets `M[a][b] = 1` for a `DirectedEdge[a, b]`. An `UndirectedEdge` also sets `M[b][a]`, so an undirected graph gives a symmetric matrix. Rows and columns follow the canonical `VertexList` order.

**Data structures.** The graph is `Graph[List[v...], List[edge...]]`. Endpoints resolve through a `GraphVIdx` hash index built once per call. This replaced a linear scan per edge that cost `O(E V)` on top of the matrix itself. The result is a freshly built `List` of `List`s of machine integers, so `Det`, `Tr` and `Eigenvalues` consume it directly. Parallel edges are forbidden by validation, so every entry is `0` or `1`.

**Complexity / limits.** `O(V^2 + E)` time and space, and the `V^2` output dominates. The matrix is dense, not sparse. A non-graph argument, or a call with more or fewer than one argument, is left unevaluated. The result is not cached.
