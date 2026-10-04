---
references:
  - "F. V. Fomin, F. Grandoni and D. Kratsch, *A measure and conquer approach for the analysis of exact algorithms*, J. ACM **56**(5) (2009) Art. 25."
source: src/graph/galg_mis.c
---
**Algorithm.** `builtin_find_independent_vertex_set` ignores edge direction. With one argument it calls the exact solver `galg_max_independent_set` and returns `{s}`, a single *maximum* independent set. With a size spec (`k`, `{k}`, `{kmin, kmax}`) and an optional count (`n` or `All`) it instead enumerates *maximal* independent sets by size through `galg_clique_spec_query`, the same engine `FindClique` uses on the complement. The solver works per connected component:

1. a bipartite component of at least 256 vertices (grids, trees, even cycles) is solved in `O(m sqrt n)` by König's theorem (`gm_konig`);
2. a dense component (average degree >= 8 or density >= 0.05, at most 3000 vertices) goes to the bitset maximum-clique branch and bound on its complement;
3. a sparse component uses branch and reduce, seeded with a greedy minimum-degree incumbent. Reductions run to a fixed point: degree 0 and 1, degree-2 triangle, degree-2 *folding* (replace `{v, a, b}` by one vertex adjacent to `N(a) u N(b)`, alpha drops by exactly 1), and domination (`N[v]` inside `N[u]` drops `u`). The bound "taken so far + greedy clique cover of the rest" prunes, a disconnected remainder splits, and otherwise the search branches on a maximum-degree vertex: drop it with its mirrors (Fomin-Grandoni-Kratsch) or take it.

**Data structures.** A `GalgUG` CSR view built from the `Graph[List, List]` expression tree (edges memoized by `graph_edge_indices`), then a `GmMis` search state: growable adjacency lists, a removal/fold trail so every reduction is undone exactly, a solution stack, and stamped scratch arrays. Dense components use bitsets. Results are rebuilt as vertex expressions from the integer positions.

**Complexity / limits.** Maximum independent set is NP-hard; the worst case is exponential. The solver gives up, and the head stays unevaluated, after 2e7 search nodes, 6e8 aggregate work, 64M ints of per-level scratch, 20000 levels of recursion, or when `TimeConstrained`'s deadline passes. It never returns a merely maximal set for the one-argument form. Which optimum is returned among ties is fixed by the vertex order and the reductions, so repeated calls agree.
