---
source: src/graph/gops_preds.c
---
**Algorithm.** `builtin_weighted_graph_q` returns `True` iff `g` is a valid graph that
carries an `EdgeWeight` option list — the single test
`graph_is_valid(g) && graph_edge_weight_list(g) != NULL`. A graph built with
`EdgeWeight -> {w1, ...}` stores that list in the canonical 3-argument form; an unweighted
graph has no such list and gives `False`. `EdgeWeightedGraphQ` is defined to call this same
function, because Mathilda has no vertex weights and so the edge-weighted and weighted
questions coincide.

**Data structures.** `graph_edge_weight_list` returns a borrowed pointer to the
`EdgeWeight` `List` stored inside the canonical `Graph[List verts, List edges, opts...]`
tree (per-edge options are read by key, not position — see `graph.h`). `graph_is_valid` is
the memoized validator, so the validity half is `O(1)` after the first query. `gops_truth`
builds the result symbol.

**Complexity / limits.** `O(1)` apart from the first validation pass. The test is purely
presence-of-option: it does not inspect whether the stored weights are numeric or
non-negative (that stricter check is `graph_weights_usable`, used by the shortest-path
code). A non-graph argument gives `False` and never stays unevaluated.
