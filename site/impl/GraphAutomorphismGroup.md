---
references:
  - "B. D. McKay and A. Piperno, *Practical graph isomorphism, II*, J. Symbolic Comput. **60** (2014) 94-112."
source: src/graph/galg_isoheads.c
---
**Algorithm.** `builtin_graph_automorphism_group` reduces the graph with `gi_build` (loops into colours, parallel edges by subdivision vertices, as in `FindGraphIsomorphism`) and calls `galg_iso_automorphisms`, which walks the individualization-refinement search tree of `galg_iso.c` with nauty-style pruning: trace pruning, automorphisms detected from leaves equal to the first or best leaf and from internal nodes whose trace matches the first path, and orbit pruning on the first path. The automorphisms found generate the full group. Each generator is turned into `Cycles[{...}]` (1-based points, each cycle starting at its least point, fixed points dropped, the identity skipped) and the result is `PermutationGroup[{gens}]`.

**Data structures.** CSR relations (`GalgIsoGraph`) built from the `Graph[List, List]` tree, a flat `int` array of generators (one permutation of the reduced vertex set per row), and a `seen` byte mask used while extracting cycles. Only generators acting on the original vertices are exposed, since subdivision vertices are determined by their ends.

**Complexity / limits.** Computing a generating set is not known to be polynomial in general, but the pruned search is fast on most graphs, and only a generating set is returned (never the elements), so the output stays compact even when the group is enormous. The search is budgeted at 5e7 nodes and respects `TimeConstrained`; on exhaustion the head stays unevaluated. A graph with no symmetry yields `PermutationGroup[{}]`. The set of generators is not unique.
