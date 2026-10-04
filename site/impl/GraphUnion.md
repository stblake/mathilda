---
source: src/graph/gops_setops.c
---
**Algorithm.** `builtin_graph_union` requires every argument to be a valid graph;
`GraphUnion[g]` returns a copy of `g`. The vertex set of the result is the union of
all the graphs' vertices in canonical (`Sort`) order. `union_build` maps each
graph's vertices to union ids — by an `O(V)` elementwise `SameQ` check when a
vertex list equals the first graph's (the common case of graphs sharing one vertex
set: no hashing at all), else through a single `GraphVIdx` hash over the union —
then orders the union with `expr_compare`, skipping the sort when it is already
sorted and running on machine integers when every vertex is one. The edges are the
**distinct** edges of all the graphs (an undirected edge equals its reversal), with
duplicates found in a `GopsKeySet` over integer edge keys. Edge order follows
Mathematica 15: all-undirected edges are oriented by canonical vertex order in
first-appearance order, all-directed edges keep first-appearance order, and a mixed
set is put in canonical order by a counting sort on `(kind, first, second)`.

**Data structures.** A `Union` struct holds the per-graph vertex-id maps and the
ordered union vertex array. Each kept edge is an `OutEdge` of result-position
endpoints plus a direction flag and the source edge node it may share. Deduplication
and canonical ordering are integer passes — a `GopsKeySet` for distinctness and a
stable counting sort (`gops_csort`) for order — so the whole assembly is `O(V + E)`
apart from the vertex sort. `union_graph` builds the result with `gops_graph_new`,
sharing a source edge node by `expr_copy` whenever its orientation is unchanged.

**Complexity / limits.** `O(V + E)` plus the `O(V log V)` vertex sort, which drops
to `O(V)` when the union is already sorted or every vertex is an integer. The result
is a simple graph; edge **weights are dropped**, as in Mathematica. `res` is
borrowed and never modified.
