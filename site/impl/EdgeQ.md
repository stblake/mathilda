---
source: src/graph/membership.c
---
**Algorithm.** `builtin_edge_q` takes `EdgeQ[g, e]`. It returns `False` when `g` is not a valid graph, or when `e` is not a two-argument edge. `query_kind` accepts `DirectedEdge` and `Rule` (`u -> v`) as directed, and `UndirectedEdge` and `TwoWayRule` (`u <-> v`) as undirected. It then calls `graph_has_edge(g, u, v, directed)`. An undirected query matches an `UndirectedEdge` in either orientation. A directed query matches only a `DirectedEdge` with the same ordered endpoints, so `1 -> 2` is not an edge of a graph whose only edge is `1 <-> 2`.

**Data structures.** The graph is `Graph[List, List]`. Membership is structural (`expr_eq`), so vertex `1` does not match `1.0`. The probe is an `O(1)` hash lookup in the validated-graph memo, which already holds the vertex index and edge-key set that validating `g` built.

**Complexity / limits.** `O(1)` per query once `g` has been validated. The first query on a new graph pays the validation cost, which is `O(V + E)`. The result is always `True` or `False`, never unevaluated, except for a wrong argument count. An edge whose endpoints are not vertices of `g` is `False`.
