---
references:
  - "S. Even, *An algorithm for determining whether the connectivity of a graph is at least k*, SIAM J. Comput. **4** (1975) 393-396."
  - "A. H. Esfahanian and S. L. Hakimi, *On computing the connectivities of graphs and digraphs*, Networks **14** (1984) 355-366."
source: src/graph/connectivity.c
---
**Algorithm.** `builtin_vertex_connectivity` gives the minimum number of vertices whose removal
disconnects the graph. The primary engine (`galg_vertex_connectivity` in `galg_flow.c`) is a
max-flow method: Even's split-vertex construction turns each vertex `v` into `v_in -> v_out`
with capacity 1, so a vertex cut becomes a minimum `s`-`t` cut, and Esfahanian–Hakimi pair
selection bounds the number of local `kappa(s,t)` computations needed for the global minimum.
Each local value is a Dinic max flow on the unit-capacity residual network. A brute-force
`C(n,k)` subset enumeration (removing every `k`-subset and testing connectivity) is retained
only as an allocation-failure fallback.

**Data structures.** The split residual network as a CSR of unit-capacity arcs, with Dinic's
BFS level array and current-arc pointers. The fallback uses a `GraphAdj`, a combination index
array, and a `removed[]` mask passed to `graph_count_components`.

**Complexity / limits.** The flow path is polynomial (each Dinic max flow on unit capacities is
about `O(E sqrt(V))`); the fallback is exponential. A complete graph returns `n-1`; an
already-disconnected graph or one with `n <= 1` returns `0`. A non-graph argument returns
unevaluated.
