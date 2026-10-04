---
source: src/graph/galg_mis.c
---
**Algorithm.** `EdgeCoverQ[g, es]` is a linear check, not a search (`gm_edge_set_q` with `cover = 1`). It gives `False` unless `g` is a valid graph and `es` is a list. Each element must be an edge of `g`: `gm_edge_of` accepts `UndirectedEdge`/`TwoWayRule` (matching in either orientation) and `DirectedEdge`/`Rule` (only as given), resolved through `graph_has_edge`. Every endpoint is marked in a visited mask, and the answer is `True` exactly when every vertex of `g` was touched. Its sibling `IndependentEdgeSetQ` shares the code and instead rejects an edge whose end was already marked.

**Data structures.** A byte array `hit[n]` indexed by vertex position, with positions looked up by `galg_vertex_arg` from the vertex list of the `Graph[List, List]` tree.

**Complexity / limits.** `O(n + |es|)` plus the edge-membership lookups. Never stays unevaluated for two arguments: a non-graph, a non-list, or an element that is not an edge of `g` gives `False`. A graph with an isolated vertex has no edge cover. Repeated edges are harmless.
