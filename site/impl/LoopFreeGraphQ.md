---
source: src/graph/gops_preds.c
---
**Algorithm.** `builtin_loop_free_graph_q` returns `gops_truth(graph_is_valid(g))`: it is
`True` for every valid Mathilda graph and `False` otherwise. A self-loop is an edge whose
two endpoints coincide, and `Graph[...]` construction rejects self-loops outright, so a
canonical graph has none — loop-freeness is guaranteed by validity, with no edge scan of
its own. (This makes `LoopFreeGraphQ` and `SimpleGraphQ` the same test in Mathilda.)

**Data structures.** The work is entirely in `graph_is_valid`, the memoized validator in
`graph_util.c`, which caches the canonical-shape / no-self-loop / no-parallel-edge result
on `g`'s node so repeat calls are `O(1)`. The graph is the ordinary
`Graph[List verts, List edges]` `Expr` tree; `gops_truth` builds the `True`/`False` symbol.

**Complexity / limits.** `O(1)` on a memo hit, one `O(V + E)` validation pass otherwise.
A non-graph argument is not valid, so the head gives `False` and never stays unevaluated.
