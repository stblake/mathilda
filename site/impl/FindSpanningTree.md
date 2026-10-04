---
references:
  - "J. B. Kruskal, *On the shortest spanning subtree of a graph and the traveling salesman problem*, Proc. Amer. Math. Soc. **7** (1956) 48-50."
  - "Y. J. Chu and T. H. Liu, *On the shortest arborescence of a directed graph*, Sci. Sinica **14** (1965) 1396-1400; J. Edmonds, *Optimum branchings*, J. Res. NBS **71B** (1967) 233-240."
source: src/graph/spanningtree.c
---
**Algorithm.** `builtin_find_spanning_tree` returns a spanning tree/forest as a
`Graph`, choosing the method from the graph's kind (checked case by case against
Mathematica 15): an **undirected unweighted** graph gives a BFS spanning forest;
an **undirected weighted** graph a minimum spanning forest by **Kruskal** (ties
broken by endpoint positions, reproducing Mathematica's choice); a **directed
unweighted** graph a BFS branching with roots taken in decreasing DFS finishing
time (fewest roots); a **directed weighted** graph a minimum-weight spanning
branching by **Chu-Liu/Edmonds** (the `O(E log V)` contraction formulation with a
rollback union-find and leftist mergeable heaps carrying a lazy additive delta).
A mixed graph is left unevaluated. `FindSpanningTree[{g, v}]` restricts the tree
to `v`'s component/reach, rooted at `v`.

**Data structures.** An out-incidence CSR with edge ids (`Inc`), and for the
weighted paths exact **GMP rational** weights (`mpq_t`): Integer/Rational/Real/
MPFR become exact binary rationals so `1/3` and `0.3333333333333333` are told
apart, other numerics go through `N[w, 40]`. Everything else runs on the memo's
integer endpoint arrays; `build_tree` re-emits edges in sorted endpoint order and
seeds the result into the graph memo.

**Complexity / limits.** `O(V + E)` unweighted, `O(E log V)` weighted. Weights
are compared exactly; a non-real weight (symbol, Complex) leaves the call
unevaluated. Trailing `Method -> ...` options are accepted and ignored — every
method yields the same optimum.
