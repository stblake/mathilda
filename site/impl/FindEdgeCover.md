---
references:
  - "T. Gallai, *Über extreme Punkt- und Kantenmengen*, Ann. Univ. Sci. Budapest. Eötvös Sect. Math. **2** (1959) 133-138."
  - "J. Edmonds, *Paths, trees, and flowers*, Canad. J. Math. **17** (1965) 449-467."
  - "J. E. Hopcroft and R. M. Karp, *An n^{5/2} algorithm for maximum matchings in bipartite graphs*, SIAM J. Comput. **2** (1973) 225-231."
source: src/graph/galg_matching.c
---
**Algorithm.** `builtin_find_edge_cover` computes a minimum edge cover — a smallest set of edges
touching every vertex — by Gallai's reduction: take a maximum matching, then add one arbitrary
incident edge for each still-uncovered vertex, giving a cover of size `n - nu(g)`. The maximum
matching starts from a Karp–Sipser greedy phase (repeatedly match a degree-1 vertex to its
neighbour) and is completed exactly by **Hopcroft–Karp** when the graph is bipartite and by
**Edmonds' blossom algorithm** otherwise (union-find blossom bases, Hungarian-tree retirement).
Edge direction is ignored, matching Mathematica's semantics.

**Data structures.** A `GalgUG` CSR, a `mate[]` matching array, Karp–Sipser degree/queue arrays,
and for the blossom path a union-find `dsu[]` with tree labels, a `dead[]` retirement mask and
LCA marking. The cover is assembled in `EdgeList` order from the original edge expressions.

**Complexity / limits.** Polynomial — `O(E sqrt(V))` on the bipartite path, and each blossom
search `O(E·alpha(V))` over the part it reaches. No hard node cap, only a `TimeConstrained`
poll. If the graph has an isolated vertex no edge cover exists and the result is `{}`. A
non-graph argument returns unevaluated. The cover returned *is* a minimum.
