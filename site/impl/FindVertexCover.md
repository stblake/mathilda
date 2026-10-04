---
references:
  - "F. V. Fomin, F. Grandoni and D. Kratsch, *A measure and conquer approach for the analysis of exact algorithms*, J. ACM **56**(5) (2009) Art. 25."
source: src/graph/galg_mis.c
---
**Algorithm.** `builtin_find_vertex_cover` takes exactly one argument and returns a *minimum* vertex cover. A set is a vertex cover exactly when its complement is independent, so it runs the same exact solver as `FindIndependentVertexSet` (`galg_max_independent_set`) and returns the vertices *not* in the maximum independent set it finds. Edge direction is ignored. The solver is per component: König's theorem for large bipartite components, the complement-clique bound for dense ones, and branch and reduce (degree-0/1/2, folding, domination, clique-cover bound, mirror branching) for sparse ones.

**Data structures.** The `Graph[List, List]` expression is flattened to a `GalgUG` CSR adjacency; a byte mask marks the independent set and the cover is read off as the unmarked vertices, in vertex order, and rebuilt as a list of the graph's own vertex expressions.

**Complexity / limits.** Minimum vertex cover is NP-hard, so the worst case is exponential. The solver gives up after a fixed node/work/memory/depth budget or when the `TimeConstrained` deadline passes, and the head then stays unevaluated; it never returns a non-minimum cover. A graph that is not valid, or any form with more than one argument, is left unevaluated. Ties between optimal covers are broken deterministically by the vertex order.
