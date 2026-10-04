---
source: src/graph/galg_mis.c
---
**Algorithm.** `VertexCoverQ[g, vlist]` is `True` when every edge of `g` has at
least one endpoint in `vlist`. It is the complement test of
`IndependentVertexSetQ` and shares the same routine `gm_vertex_set_q` (called
with `cover = 1`): `gm_vertex_positions` resolves the list to vertex indices
(`False` if any element is not a vertex of `g`), the chosen vertices are marked in
a flag array, and one pass over the edge-index arrays `(eu, ev)` reports `False`
as soon as an edge has **neither** endpoint in the set. Edge direction is ignored
and repeated vertices are allowed, matching Wolfram's semantics. (A vertex cover
and a maximum independent set are exact complements, which is why
`FindVertexCover` in the same file delegates to the independent-set solver and
returns the complement.)

**Data structures.** Endpoints are read from the graph's cached index arrays via
`graph_edge_indices`; membership is a `char in[n]` bitmap and the resolved
positions a small `int` array. No `GalgUG` or bitset is built — a single linear
edge scan suffices — and all scratch is freed on every path.

**Complexity / limits.** `O(|vlist| + m)` time, `O(n)` space. The predicate is
total: it returns a Boolean for a non-graph argument or a non-vertex list element
(both `False`) and never stays unevaluated.
