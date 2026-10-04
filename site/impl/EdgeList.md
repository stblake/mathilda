---
source: src/graph/edgelist.c
---
**Algorithm.** `builtin_edge_list` validates `g` and returns a deep copy of the
graph's canonical edge sublist `args[1]`, preserving stored order. Each element
is already in the canonical two-argument `DirectedEdge[u, v]` / `UndirectedEdge[u, v]`
form the constructor normalises edge sugar (`u -> v`, `u <-> v`) into, so
`EdgeList` is the inverse of the edge half of `Graph[...]`. When the argument is
not a valid graph the call is forwarded to `hyp_edge_list` (which returns the
hyperedges of a `Hypergraph`) and otherwise yields `NULL`.

**Data structures.** The result is `expr_copy(g->data.function.args[1])` — an
independent `List` of edge nodes, so the caller can mutate or consume it without
disturbing the graph. No adjacency structure is consulted.

**Complexity / limits.** `O(E)` for the copy, after an `O(1)` memo-hit (or
`O(V + E)` first-time) validation. One-argument only; other arities return
`NULL`.
