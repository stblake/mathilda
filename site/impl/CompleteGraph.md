---
source: src/graph/generators.c
---
**Algorithm.** `builtin_complete_graph` builds the undirected complete graph
`K_n` on the vertices `1..n`: it emits every pair `i < j` as an
`UndirectedEdge[i, j]`, giving exactly `n(n-1)/2` edges. The multipartite form
`CompleteGraph[{n1, n2, ...}]` is handled by a second registration,
`builtin_gmet_complete_graph` in `src/graph/gmet_generators.c`, which wraps the
base builtin: it lays the parts out as consecutive vertex blocks and joins every
pair of vertices in *different* parts, so `CompleteGraph[{2, 3}]` is the complete
bipartite graph on `2 + 3` vertices.

**Data structures.** Both forms assemble plain C arrays of `Expr*` for the
vertices and edges and hand them to `expr_new_function` as a
`Graph[List, List]`; the evaluator's `builtin_graph` canonicalizes and validates
the result and seeds the graph memo. Edges are generated in row-major pair order
(`undirected_edge(i, j)` shares a single interned `UndirectedEdge` head path).

**Complexity / limits.** `O(n^2)` — the output has `Theta(n^2)` edges, so the
cost is output-bound. A negative or non-integer `n` leaves the call unevaluated.
