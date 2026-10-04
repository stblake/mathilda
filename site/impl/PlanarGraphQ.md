---
references:
  - "U. Brandes, *The Left-Right Planarity Test*, manuscript (2009), algorithms 2-5."
source: src/graph/galg_planar.c
---
**Algorithm.** `builtin_planar_graph_q` gives `False` for a non-graph and otherwise runs `galg_planar_test` on the underlying simple undirected graph (direction and anti-parallel pairs do not matter). It is the left-right planarity test of de Fraysseix and Rosenstiehl in Brandes's formulation. Trivial cases first: at most 4 vertices or fewer than 9 edges is planar (K5 and K3,3 need 9), and `m > 3n - 6` is rejected by the Euler bound. Then per connected component two DFS passes: *orientation* orients edges away from the root and computes `lowpt`, `lowpt2` and a nesting depth `2*lowpt + [lowpt2 < height]`, ordering each vertex's outgoing edges by it; *testing* walks the tree in nesting order keeping a stack of conflict pairs (L, R) of return-edge intervals, merging constraints after each child (`add_constraints`) and trimming returns when retreating (`remove_back_edges`). A contradiction means non-planar. The embedding phase is omitted.

**Data structures.** CSR adjacency from the `Graph[List, List]` tree; one malloc'd block of ints (heights, parent edges, lowpoints, `ref` chains, a flat conflict-pair stack of four edge ids per pair). Both DFS passes are iterative with an explicit stack and per-vertex cursors, so a path of 10^6 vertices needs no call stack; phase 2 runs on preorder numbers for cache locality.

**Complexity / limits.** `O(n + m)` time and memory. Returns a definite `True`/`False`; only an allocation failure leaves the head unevaluated. The test decides planarity but constructs no embedding or Kuratowski subgraph.
