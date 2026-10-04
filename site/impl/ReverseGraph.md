---
source: src/graph/gops_transform.c
---
**Algorithm.** `builtin_reverse_graph` reverses every directed edge of `g` while
leaving undirected edges, the edge order, and the weights untouched. When `g` has
no directed edges it returns a copy of the input unchanged (a fast exit). For a
graph with directed edges it walks the edge list and, for each `DirectedEdge[a,
b]`, emits `DirectedEdge[b, a]` with its endpoints swapped; an undirected edge is
copied verbatim. The `k`-th output edge corresponds to the `k`-th input edge, so
`EdgeList` order is preserved and the `EdgeWeight`/`EdgeCapacity` lists stay
aligned.

**Data structures.** A `GopsView` of the input and fresh vertex, edge, weight,
and integer-endpoint arrays; reversed edges share a single interned
`DirectedEdge` head (`GopsHeads`). `gops_graph_new` seeds the result into the
graph memo.

**Complexity / limits.** `O(V + E)`, a single linear pass. Weights are carried
over unchanged, so a weighted reversed graph keeps its weights. A non-graph
argument leaves the call unevaluated.
