---
source: src/graph/graphprops.c
---
**Algorithm.** `builtin_bipartite_graph_q` tests whether the *underlying
undirected* graph is 2-colourable. It runs a breadth-first 2-colouring over
every component: each vertex is given a side (0/1) and a neighbour on the same
side is a conflict, which proves an odd cycle and so non-bipartiteness. Edge
direction is ignored — the neighbourhood of a vertex is both its successors
(`out[]`) and its predecessors (`in[]`) — because bipartiteness is a property of
the vertex split, not of orientation. An edgeless graph is trivially bipartite;
a non-graph argument gives `False` (never unevaluated), as every `*Q` predicate
does.

**Data structures.** It builds the integer-indexed `GraphAdj` (CSR successor and
predecessor lists, `src/graph/graph.h`), a `signed char side[]` array, and an
explicit BFS queue `q[]`. The result is memoized on the graph node through
`graph_prop_set(g, GRAPH_PROP_BIPARTITE, ...)`, so a repeat query on the same
`Graph` object is `O(1)` — the same property-caching an atomic Mathematica
`Graph` does.

**Complexity / limits.** `O(V + E)` for the first call, `O(1)` thereafter.
Allocation failure is the only path that returns `NULL`; a malformed graph
returns `False`.
