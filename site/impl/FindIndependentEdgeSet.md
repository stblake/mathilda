---
references:
  - "J. E. Hopcroft and R. M. Karp, *An n^{5/2} algorithm for maximum matchings in bipartite graphs*, SIAM J. Comput. **2** (1973) 225-231."
  - "J. Edmonds, *Paths, trees, and flowers*, Canad. J. Math. **17** (1965) 449-467."
  - "R. M. Karp and M. Sipser, *Maximum matchings in sparse random graphs*, FOCS 1981, 364-375."
source: src/graph/galg_matching.c
---
**Algorithm.** `FindIndependentEdgeSet[g]` returns a maximum **matching** (a
largest set of pairwise non-adjacent edges), ignoring edge direction as Wolfram
does. `galg_max_matching` starts from a Karp-Sipser greedy matching
(`gm_karp_sipser`): it repeatedly matches a degree-1 vertex to its only neighbour
— always safe, some maximum matching contains that edge — and otherwise an
arbitrary remaining edge, which on sparse graphs is already optimal or within a
few edges. It then completes to an exact maximum: a bipartiteness BFS
(`gm_bipartite`) routes bipartite graphs to Hopcroft-Karp (`gm_hopcroft_karp`,
phases of shortest vertex-disjoint augmenting paths) and everything else to
Edmonds' blossom algorithm (`gm_edmonds`), one alternating-tree search per
exposed vertex with odd cycles contracted into blossoms.

**Data structures.** The graph is a CSR `GalgUG`; the matching is a single
`mate[]` array (`mate[v]` = partner or `-1`). Blossom contraction keeps bases in
a union-find forest (`gm_find` with path halving), so one search costs
`O(E·α(V))` over the part of the graph it reaches; only the vertices a search
touched are reset afterwards, and a search that fails leaves a Hungarian tree
whose vertices are retired for the rest of the run (Edmonds' lemma). The matched
edges are emitted as `g`'s own edge expressions in `EdgeList` order by
`gm_matched_edges`.

**Complexity / limits.** Hopcroft-Karp is `O(E√V)`; the blossom search is
`O(V·E·α)` overall. Both are exact maximum-cardinality algorithms, so the result
size is the matching number `ν(g)`; the particular maximum matching is not
specified (the greedy start and search order fix it). The companion
`FindEdgeCover` reuses the matching through the Gallai identity (a maximum
matching plus one incident edge per exposed vertex, size `n − ν(g)`).
