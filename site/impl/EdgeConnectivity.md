---
references:
  - "H. Nagamochi and T. Ibaraki, *Computing edge-connectivity in multigraphs and capacitated graphs*, SIAM J. Discrete Math. **5** (1992) 54-66."
  - "E. A. Dinic, *Algorithm for solution of a problem of maximum flow in a network with power estimation*, Soviet Math. Dokl. **11** (1970) 1277-1280."
source: src/graph/galg_flow.c
---
**Algorithm.** `builtin_edge_connectivity` gives the minimum total capacity of a
set of edges whose removal disconnects `g` (`EdgeWeight` as the capacity when
present, else `1` per edge). `EdgeConnectivity[g]` is the *global* minimum cut
and `EdgeConnectivity[g, s, t]` the minimum `s-t` cut. For an undirected graph
the global cut uses **Nagamochi-Ibaraki**: each round computes a
maximum-adjacency order, lowers the bound `lambda` to the least weighted degree,
and contracts every edge whose MA-attachment reaches `lambda` (so it cannot lie
on a lighter cut) — it returns Stoer-Wagner's value while usually contracting
most of the graph per round. A directed global cut is the minimum over `v != v0`
of the flows `v0 -> v` and `v -> v0`. The `s-t` form is a single max flow.

**Data structures.** Capacities are scaled to exact `int64` in a `GfCap` (integer
capacities stay exact; reals scale by a common power of two; `Infinity` gets a
finite stand-in above every finite sum). Flows run Dinic on a CSR residual
network (`GfNet`: `first[]`/`to[]`/`rev[]`/`res[]` arcs, reverse arcs adjacent),
with level-BFS and current-arc blocking-flow search. The undirected cut keeps its
own CSR plus a lazy max-heap for the MA order.

**Complexity / limits.** All arithmetic is on `int64`; the result is an `Integer`
for integer capacities and a `Real` otherwise. A graph with `< 2` vertices is
left unevaluated, as is a negative or symbolic capacity.
