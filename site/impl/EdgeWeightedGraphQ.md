---
source: src/graph/gops_preds.c
---
**Algorithm.** `builtin_edge_weighted_graph_q` is defined as a direct call to
`builtin_weighted_graph_q`, so it answers `True` iff `g` is a valid graph carrying an
`EdgeWeight` option list (`graph_is_valid(g) && graph_edge_weight_list(g) != NULL`). The
two predicates are deliberately identical: Mathilda has no vertex weights, so "edge
weighted" and "weighted" are the same property.

**Data structures.** Delegation means there is no separate state — `graph_edge_weight_list`
returns a borrowed pointer to the `EdgeWeight` `List` inside the canonical
`Graph[List verts, List edges, opts...]` tree, and `graph_is_valid` is the memoized
validator, so the validity half is `O(1)` after the first query on a node.

**Complexity / limits.** `O(1)` apart from the first validation pass. Like `WeightedGraphQ`
it tests only that an `EdgeWeight` list is present, not that the stored weights are numeric
or non-negative. A non-graph argument gives `False` and never stays unevaluated.
