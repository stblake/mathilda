---
source: src/graph/generators.c
---
**Algorithm.** `builtin_cycle_graph` builds the undirected cycle `C_n` on the
integer vertices `1..n`. It emits the path edges `1 <-> 2, 2 <-> 3, ...,
(n-1) <-> n`, then adds the wrap edge `n <-> 1` only when `n >= 3`: for `n <= 2`
that edge would duplicate one already present (`n = 2` is a single edge, `n <= 1`
is edgeless), and Mathilda graphs are simple. The assembled
`Graph[List verts, List edges]` is returned and the evaluator canonicalises and
validates it through `builtin_graph`. A non-integer or negative argument, or any
arity other than one, returns `NULL` and the call is left unevaluated.

**Data structures.** Vertices are `EXPR_INTEGER` nodes `1..n` from `int_vertices`;
the edges are `UndirectedEdge` nodes in a `calloc`'d `Expr*` array. `make_graph`
wraps both arrays into the two `List`s of a `Graph[...]`, moving ownership of the
array contents into the new nodes (the arrays themselves are then freed). No
adjacency or incidence structure is built here — that is the validated-graph memo's
job, populated lazily once `builtin_graph` accepts the tree.

**Complexity / limits.** `O(n)` vertices and `O(n)` edges in a single allocation
pass. The output is a simple undirected graph — the `n <= 2` wrap-edge guard is
exactly what keeps it so — in which every vertex has degree two for `n >= 3`. `res`
is freed by the evaluator on success.
