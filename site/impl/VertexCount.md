---
source: src/graph/counts.c
---
**Algorithm.** `builtin_vertex_count` returns the number of vertices of `g` as an
integer — a thin reader over the canonical form `Graph[List verts, List edges]`:
the answer is the argument count of the vertex list. When the argument is not a
valid graph it falls through to `hyp_vertex_count`, which handles a `Hypergraph`
(and otherwise returns `NULL`, leaving the expression unevaluated). The companion
`EdgeCount` is the identical reader over the edge list.

**Data structures.** None — it reads `g->data.function.args[0]->...arg_count`
after `graph_is_valid` confirms the shape and populates the memo. The result is a
fresh `Integer`.

**Complexity / limits.** `O(1)` for a canonical graph (plus the one-time
`O(V + E)` validation on first contact with the node). Returns `NULL` for a
non-graph, non-hypergraph argument.
