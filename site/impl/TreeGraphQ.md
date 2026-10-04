---
source: src/graph/acyclic.c
---
**Algorithm.** `builtin_tree_graph_q` returns `False` for a non-graph, then reads the cached `GRAPH_PROP_TREE` flag. On a miss it tests that the graph has at least one vertex and exactly `V - 1` edges, then runs union-find over all edges. It answers `True` when every edge joins two different components. With `V - 1` joins the graph is one component, so it is connected and acyclic. Direction is ignored, so the out-tree `1 -> 2, 1 -> 3` is a tree. The anti-parallel pair `1 -> 2, 2 -> 1` is two edges on two vertices and is not. The null graph is not a tree.

**Data structures.** It reads the memoized `graph_edge_indices` views (`eu[k]`, `ev[k]`, `directed[k]`) over the `Graph[List, List]` expression. Union-find uses two `int` arrays (`parent` and `size`) with path halving and union by size. The answer is stored with `graph_prop_set`.

**Complexity / limits.** `O(V + E alpha(V))` on the first query, and `O(1)` on a repeat for the same graph node. The edge count test rejects most non-trees before any union-find work. The result is always `True` or `False`: an argument that is not a valid graph gives `False`, not an unevaluated form.
