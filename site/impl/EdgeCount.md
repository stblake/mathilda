---
source: src/graph/counts.c
---
**Algorithm.** `builtin_edge_count` is a thin reader over the canonical form
`Graph[List[verts], List[edges]]`: it validates `g` and returns the integer
`arg_count` of the edge sublist `args[1]`. Every edge in the canonical list —
whether a `DirectedEdge` or an `UndirectedEdge` — is counted exactly once, so
the count is the length of `EdgeList[g]`. When the argument is not a valid graph
the call is forwarded to `hyp_edge_count`, which answers for a `Hypergraph` and
otherwise returns `NULL` (the call stays unevaluated).

**Data structures.** None beyond the result. The edge list is already a stored
child of the graph node, so no traversal or allocation happens past constructing
the fresh integer the evaluator then owns.

**Complexity / limits.** `O(1)` — the arity is read straight from the stored
list header. The one-time validation `graph_is_valid(g)` is `O(1)` on a memo hit
and `O(V + E)` the first time a graph is seen. This head takes exactly one
argument; any other arity returns `NULL`, so there is no edge-pattern counting
form (`EdgeCount[g, patt]` is left unevaluated).
