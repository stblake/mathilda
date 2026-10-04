---
source: src/graph/gops_preds.c
---
**Algorithm.** `builtin_simple_graph_q` returns `gops_truth(graph_is_valid(g))` — it is
`True` for every valid Mathilda graph and `False` otherwise. The reason is structural:
`Graph[...]` construction rejects self-loops and parallel/duplicate edges, so a canonical
graph is *by construction* simple. There is no separate scan for multi-edges or loops to
run; validity is the whole question.

**Data structures.** `graph_is_valid` is the memoized validator from `graph_util.c`: on
`g`'s first use it checks the canonical shape (`Graph[List verts, List edges]`, every edge
a 2-argument `DirectedEdge`/`UndirectedEdge`, no self-loops, no parallel edges, endpoints
all in `verts`) and caches the result on the node, so repeat calls are `O(1)`. `gops_truth`
maps the `int` to a fresh `True`/`False` symbol.

**Complexity / limits.** `O(1)` on a memo hit, one `O(V + E)` validation pass otherwise.
A non-graph argument is simply not valid, so the head gives `False`; it never stays
unevaluated.
