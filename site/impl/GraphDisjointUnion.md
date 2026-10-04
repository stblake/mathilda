---
source: src/graph/gops_setops.c
---
**Algorithm.** `builtin_graph_disjoint_union` places the input graphs side by
side on fresh integer vertices `1..n`: graph `g1`'s vertices first (in
`VertexList` order), then `g2`'s, and so on, with each graph's edges translated
by its block offset. No vertices are identified even when the inputs share labels
— this is exactly what distinguishes it from `GraphUnion`, which merges equal
vertices. `GraphDisjointUnion[g]` is `g`. Every edge keeps its direction; weights
are dropped.

**Data structures.** It first sums the vertex and edge counts to size the output
arrays, fills the relabelled vertex list `{1, ..., n}`, then opens a `GopsView`
on each input and copies its edges with endpoints shifted by a running `off`.
Edge nodes are built through a shared interned-head cache (`GopsHeads`) and the
result goes through `gops_graph_new`, which seeds the graph memo.

**Complexity / limits.** `O(V + E)` over all inputs — a single linear pass, no
hashing (unlike the vertex-union set operations). Each argument must be a valid
graph, else the call is unevaluated.
