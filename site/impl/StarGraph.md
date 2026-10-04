---
source: src/graph/generators.c
---
**Algorithm.** `builtin_star_graph[n]` reads a non-negative integer count via `as_count` and joins hub vertex 1 to each of vertices `2..n`, giving exactly `n - 1` edges `UndirectedEdge[1, i]`. `StarGraph[1]` is one isolated vertex and `StarGraph[0]` is the empty graph.

**Data structures.** A `calloc`'d array of `n - 1` edge nodes built by `undirected_edge` and an integer vertex list from `int_vertices`, passed to `make_graph`, which assembles the canonical `Graph[List[1..n], List[edges]]` `Expr` tree. No duplicate edge is possible, so no dedup pass is needed (unlike `CycleGraph`'s wrap edge).

**Complexity / limits.** `O(n)` time and memory. A non-integer or negative argument, or a wrong argument count, leaves the call unevaluated.
