---
source: src/graph/gops_preds.c
---
**Algorithm.** `builtin_mixed_graph_q` is `True` exactly when `g` has **both** a
directed and an undirected edge. It reads the directed-edge count from the memo
(`graph_directed_edge_count`) and compares it with the total edge count: the
graph is mixed iff the count is strictly between `0` and the total. A graph with
no directed edges (all undirected or edgeless) or with every edge directed is not
mixed.

**Data structures.** None of its own — two integers off the validated-graph memo
entry. The result is a fresh `True`/`False` symbol.

**Complexity / limits.** `O(1)` on a memoized graph. A non-graph argument gives
`False`, never unevaluated — as every `*Q` predicate does. Several transforms
(`UndirectedGraph`, `DirectedGraph`, `LineGraph`, `FindSpanningTree`) decline on
mixed graphs, so this is the predicate that tells them apart.
