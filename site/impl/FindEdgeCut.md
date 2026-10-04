---
references:
  - "H. Nagamochi and T. Ibaraki, *Computing edge-connectivity in multigraphs and capacitated graphs*, SIAM J. Discrete Math. **5** (1992) 54-66."
  - "E. A. Dinic, *Algorithm for solution of a problem of maximum flow in a network with power estimation*, Soviet Math. Dokl. **11** (1970) 1277-1280."
source: src/graph/galg_flow.c
---
**Algorithm.** `builtin_find_edge_cut` returns the *edges* of a minimum cut (in
`EdgeList` order), rather than its weight. `FindEdgeCut[g]` finds a global
minimum cut — Nagamochi-Ibaraki for an undirected graph, the min over `s-t`
flows for a directed one — then reports the edges crossing the shore boundary
(`gf_cut_edges`: an undirected edge whenever its endpoints lie on opposite
sides, a directed one only when it runs side-1 to side-0). `FindEdgeCut[g, s, t]`
runs a single Dinic max flow from `s` to `t` and takes the cut closest to `s`
(the set of vertices still residual-reachable from `s`). `EdgeWeight` is used as
the capacity when present, else `1`.

**Data structures.** The same `GfCap` scaled-capacity parse and `GfNet` CSR
residual network as `EdgeConnectivity` and `FindMaximumFlow`; the shore is a
`char side[]` filled by a residual-reachability BFS (`gf_reach_from`), from which
the crossing edges are collected over the memo's integer endpoint arrays.

**Complexity / limits.** Dominated by the flow/contraction phase. The `s-t` form
returns `{}` when `s` and `t` are already separated by an empty cut; `FindEdgeCut`
of a graph with one vertex is `{}`. A negative or symbolic capacity leaves the
call unevaluated.
