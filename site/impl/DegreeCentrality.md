---
source: src/graph/gmet_centrality.c
---
**Algorithm.** `builtin_degree_centrality` accepts `[g]`, `[g, "In"]` or `[g, "Out"]`. One pass over the edge list increments `out[u]` and `in[v]` for each edge, and for an undirected edge also `out[v]` and `in[u]`, so an undirected edge counts as one arc each way. `"In"` and `"Out"` return the respective counts. The default returns their sum, except on a purely undirected graph, where `out` already equals the degree. In a mixed graph an undirected edge therefore adds 2 to the total at each endpoint, as in Mathematica. Weights are ignored and the result is exact integers.

**Data structures.** The edge-index view from `graph_edge_indices` supplies the parallel `eu`, `ev` and `directed` arrays. Two `int64_t` count arrays of length `n` are the only working storage, and the answer is a packed integer vector (`gmet_int_vector`).

**Complexity / limits.** `O(n + m)` time and `O(n)` space. Any mode string other than `"In"` or `"Out"`, or a non-string second argument, leaves the call unevaluated.
