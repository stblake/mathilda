---
references:
  - "H. Nagamochi and T. Ibaraki, *Computing edge-connectivity in multigraphs and capacitated graphs*, SIAM J. Discrete Math. **5** (1992) 54-66."
  - "M. Stoer and F. Wagner, *A simple min-cut algorithm*, J. ACM **44** (1997) 585-591."
  - "E. A. Dinic, *Algorithm for solution of a problem of maximum flow in networks with power estimation*, Soviet Math. Dokl. **11** (1970) 1277-1280."
source: src/graph/galg_flow.c
---
**Algorithm.** `builtin_find_minimum_cut` computes a **global** minimum edge cut
of `g` and returns `{value, {part1, part2}}`. Capacities are the `EdgeWeight`
when `g` carries one, else `1`. An all-undirected graph goes to
`gf_ni_mincut`, the Nagamochi-Ibaraki round structure: each round builds a
maximum-adjacency order, lowers the global bound `λ` to the least weighted
degree seen, and contracts every edge whose MA attachment `q(e) ≥ λ` (such edges
carry `λ(u,v) ≥ λ`, so contracting them cannot destroy a lighter cut) together
with the last two vertices of the order. It yields exactly Stoer-Wagner's answer
but typically collapses most of the graph per round. A graph with any directed
edge goes to `gf_directed_mincut`: the minimum over every `v ≠ v₀` of the two
`s-t` max flows `v₀→v` and `v→v₀`, each computed by Dinic's blocking-flow
algorithm (`gf_dinic`) and capped at the running bound.

**Data structures.** Non-integer capacities are scaled by a common power of two
so all arithmetic stays exact on `int64` (a machine real becomes an exact dyadic
integer; a rational goes through its `double`, the precision Mathematica reports).
The flow path uses a CSR residual network (`GfNet`) whose forward and reverse
arcs of one edge are adjacent pairs, with BFS levels and current-arc pointers;
the final cut shore is recovered by residual reachability from the source. The
builtin collects the two shores into `int` arrays and renders them with
`galg_vertex_list` in `VertexList` order.

**Complexity / limits.** Max flow is Dinic (`O(E√E)` on unit capacities, a handful
of phases in practice); the directed global cut runs `O(n)` such flows, the
undirected one is the near-linear-per-round Nagamochi-Ibaraki contraction.
Integer capacities give an exact `Integer` value, rational or real ones a `Real`,
and `Infinity` is an allowed capacity. A graph with fewer than two vertices, or a
negative or symbolic capacity, leaves the call unevaluated; the source side is
listed first for a directed graph, the side without `VertexList[g]⟦1⟧` first for
an undirected one.
