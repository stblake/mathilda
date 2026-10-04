---
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hypergraph_edge_add` appends hyperedges to `h`. `edit_items`
with `edges = 1` reads the second argument as a *list of hyperedges* when it is a
List whose every element is itself a List, otherwise as a single hyperedge. Each
appended hyperedge's vertices that are new (not in the memo's vertex index, and
de-duplicated among the additions via a temporary `GraphVIdx extra`) are added to
the vertex set in first-appearance order; the hyperedges are then appended after
the existing ones. Repeats are allowed — the result is a multi-hypergraph.

**Data structures.** The memo's vertex `GraphVIdx` plus a temporary `extra` index;
the vertex and hyperedge Lists are rebuilt with `mk_hyp`/`mk_list` (existing
hyperedges shared by `expr_copy`, appended ones copied from the argument).

**Complexity / limits.** `O(n + Σ|existing| + Σ|added|)`. A plain non-List element
inside a would-be hyperedge leaves the call unevaluated.
