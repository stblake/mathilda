---
source: src/graph/gops_edit.c
---
**Algorithm.** `builtin_subgraph` gives the subgraph of `g` **induced** by a set
of vertices: every edge of `g` with both endpoints in the set is kept. The second
argument may be an explicit vertex list (elements that are not vertices of `g`
are ignored; an edge element makes it decline, since edge-induced subgraphs are
not handled), a pattern (every matching vertex is selected via `gops_matchq`), or
a single vertex. Selected vertices are given result positions in the order
listed, then `build_induced` keeps each edge whose endpoints both survive,
emitting edges in lower-triangular adjacency order of that vertex order. Edge
weights are kept.

**Data structures.** A `pos[]` map from original to result vertex index, a
`GopsView` of the graph, and an induced-subgraph builder that walks each kept
vertex's incidence, carrying the aligned `EdgeWeight` entries and sharing the
edge nodes. The result is seeded into the memo (`gops_graph_new`).

**Complexity / limits.** `O(V + E)` plus one position lookup per selected vertex.
A single-vertex argument that is not a vertex of `g`, or an edge-list argument,
leaves the call unevaluated. The companion `NeighborhoodGraph` induces on a
distance ball instead of an explicit set.
