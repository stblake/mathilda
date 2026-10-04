---
source: src/graph/edgeweight.c
---
**Algorithm.** `builtin_edge_weight` returns `g`'s per-edge weights in `EdgeList`
order. It delegates to the shared `graph_resolve_edge_weights`
(`src/graph/graph_util.c`): when `g` carries an `EdgeWeight -> {w1, ...}` option
(the 3-argument canonical form) it returns a copy of that list, and otherwise
`{1, 1, ..., 1}` — one `1` per edge. Treating an unweighted edge as weight `1`
matches the Wolfram Language and means every valid graph has a well-defined
weight list, not just those built with explicit weights.

**Data structures.** The weight list lives inside the canonical `Graph` as a
borrowed `EdgeWeight` option list (read via `graph_edge_weight_list`, never by
argument position — which option is present decides where it sits, per
`src/graph/graph.h`). The returned list is a fresh copy so the caller owns it.

**Complexity / limits.** `O(E)`. Sharing the resolver with
`WeightedAdjacencyMatrix` guarantees the two builtins never disagree on what
"unweighted" defaults to. Returns `NULL` (unevaluated) when `g` is not a valid
graph.
