---
source: src/graph/gops_transform.c
---
**Algorithm.** `builtin_undirected_graph` gives the undirected graph underlying
`g`: each edge loses its direction and the opposite arcs `u -> v` and `v -> u`
merge into a single `u <-> v` whose weight is the **sum** of theirs. An undirected
`g` is returned unchanged. The builder maps each edge to a normalized `(lo, hi)`
endpoint pair (lower `VertexList` position first), stably counting-sorts the
edges by `(lo, hi)` (row-major, upper triangle — Mathematica's order), then walks
runs of equal pairs: a run of one keeps its weight, a longer run has its weights
summed by building and evaluating a `Plus`. Because summing weights can run
arbitrary code and evict `g` from the memo, it first takes private endpoint
copies (`gops_view_own`).

**Data structures.** A `GopsView` (owned when weighted), per-edge `perm`/`lo`/`hi`
arrays and the subsystem counting sort (`gops_csort`), fresh vertex/edge/weight/
endpoint arrays, and a shared interned `UndirectedEdge` head. `gops_graph_new`
seeds the result.

**Complexity / limits.** `O(V + E)` apart from the stable sort. The opposite
form, `DirectedGraph`, replaces each `u <-> v` by the arc pair `u -> v`, `v ->
u`.
