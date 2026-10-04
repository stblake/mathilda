---
source: src/graph/galg_mis.c
---
**Algorithm.** `IndependentVertexSetQ[g, vlist]` is `True` when no two vertices of
`vlist` are adjacent in `g`. It and `VertexCoverQ` share one routine,
`gm_vertex_set_q` (called with `cover = 0`): `gm_vertex_positions` resolves each
element of `vlist` to a vertex index (returning `False` if any element is not a
vertex of `g`), a membership flag array `in[]` marks the chosen vertices, and a
single pass over the edge-index arrays `(eu, ev)` reports `False` as soon as an
edge has **both** endpoints in the set. Edge direction is ignored (Wolfram's
semantics), so an `UndirectedEdge` and a `DirectedEdge` both count as an
adjacency, and repeated vertices in `vlist` are allowed.

**Data structures.** The edge endpoints come straight from the graph's cached
index arrays via `graph_edge_indices` — no `GalgUG` or bitset is built, since the
predicate only needs one linear edge scan. Membership is a `char in[n]` bitmap
and the chosen positions are a small `int` array; both are freed on every path.

**Complexity / limits.** `O(|vlist| + m)` time and `O(n)` space: one scan to mark
the set, one scan over the edges. The predicate is total — it returns a Boolean
for a non-graph first argument or a non-vertex list element (both `False`) and
never stays unevaluated, as the `*Q` convention requires.
