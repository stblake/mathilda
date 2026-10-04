---
references:
  - "S. Even, *An algorithm for determining whether the connectivity of a graph is at least k*, SIAM J. Comput. **4** (1975) 393-396."
  - "E. A. Dinic, *Algorithm for solution of a problem of maximum flow in networks with power estimation*, Soviet Math. Dokl. **11** (1970) 1277-1280."
source: src/graph/galg_flow.c
---
**Algorithm.** `builtin_find_vertex_cut` returns an actual minimum set of vertices whose removal
disconnects the underlying undirected graph (or, in the three-argument form, separates `s` from
`t`). It uses the same Even split-vertex construction and Esfahanian–Hakimi pair selection as
`VertexConnectivity`, but then extracts the separator from the Dinic max flow: the cut vertices
are those whose in-copy can no longer reach `t_in` while their out-copy still can. Ties are
broken toward the cut closest to `t`.

**Data structures.** A sorted-CSR undirected view `GalgUG` of the graph, the split residual
network `GfNet`, a `sep[]` output index array, and a `mark[]` reachability array over the
residual graph after the final flow.

**Complexity / limits.** Polynomial — Dinic max flow on unit capacities, as for
`VertexConnectivity`. `FindVertexCut[g]` on a complete undirected graph returns its first `n-1`
vertices; `FindVertexCut[g, s, t]` with `s` and `t` adjacent returns `{}` (they cannot be
separated without deleting an endpoint). A non-graph argument, or `s == t`, returns unevaluated.
