---
source: src/graph/connectivity.c
---
**Algorithm.** `builtin_connected_graph_q` follows Mathematica's rule: a graph
carrying any directed edge is "connected" only when **strongly** connected, and a
purely undirected graph when it has a single component. It builds the integer
adjacency and branches: with at least one directed edge it labels strongly
connected components (`graph_strong_label`, iterative Tarjan over the
out-adjacency, where an undirected edge links both ways) and requires exactly one
component; otherwise it counts undirected components (`graph_count_components`,
union/BFS). The empty graph is not connected (`n >= 1` is required).

**Data structures.** The `GraphAdj` CSR successor/predecessor lists
(`src/graph/graph.h`) plus a `comp[]` labelling array for the strong-component
pass. Everything runs on integer vertex indices from the validated-graph memo;
no expression is hashed during the scan.

**Complexity / limits.** `O(V + E)`. `NULL` (unevaluated) only on allocation
failure; a non-graph argument also returns `NULL`, since the builtin requires a
buildable adjacency.
