---
source: src/graph/adjgraph.c
---
**Algorithm.** `builtin_adjacency_graph` is the inverse of `AdjacencyMatrix`: it
reads an `n x n` matrix of literal `0`/`1` integers and builds a graph on the
vertices `1..n`. It first tests the matrix for symmetry; a symmetric matrix
yields an **undirected** graph with one `UndirectedEdge[i, j]` per pair `i < j`
with `m[i][j] = 1`, and an asymmetric matrix yields a **directed** graph with a
`DirectedEdge[i, j]` for each off-diagonal `m[i][j] = 1`. Diagonal entries
(self-loops) are ignored. Any entry that is not `0` or `1`, or a non-square or
non-list argument, leaves the call unevaluated (returns `NULL`).

**Data structures.** The matrix rows are read directly as `Expr` lists. The
builder allocates the vertex list `{1, ..., n}` and an upper-bounded edge array
(`n*n` slots, trimmed by the actual count), then wraps them in a
`Graph[List, List]` expression. That raw form is handed back to the evaluator,
whose `builtin_graph` normalizes and validates it and seeds the per-node graph
memo (the shared `GraphVIdx`/endpoint arrays described in `src/graph/graph.h`).

**Complexity / limits.** `O(n^2)` — the symmetry test and the fill are both full
passes over the matrix. The round-trip `AdjacencyGraph[AdjacencyMatrix[g]]`
reproduces `g` exactly when `g`'s vertices are `1..n`. Weights are not read: the
input is a plain 0/1 adjacency matrix, so the result is unweighted.
