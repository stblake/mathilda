---
source: src/graph/gops_edit.c
---
**Algorithm.** `builtin_edge_delete` removes matching edges and rebuilds the graph on the *same*
vertex set (no vertices are dropped). Each argument item is classified individually: an item
containing a pattern is matched against every edge and clears the keep-flag of each hit; a
literal edge is resolved to its `EdgeList` index (and the call stays unevaluated if that edge is
not in the graph). An undirected edge matches either orientation, because edge keys canonicalise
undirected endpoints. The surviving edges and all vertices are copied into a fresh canonical
`Graph`, remapping endpoint indices.

**Data structures.** A `GopsView` exposes the raw vertex/edge arrays plus integer endpoints, so
the literal-edge pass is an integer pass against a `GopsKeySet` (an open-addressing hash of the
argument's edges, probed once per graph edge). The view is taken *after* any pattern matching,
since matching can evaluate arbitrary user code. `EdgeWeight` entries are kept aligned with the
surviving edges.

**Complexity / limits.** `O(V + E)` plus one hash probe per edge for the literal path; the
pattern path costs `O(items·E)` matches. The result is a canonical `Graph`; a non-graph argument,
or a literal that is not an edge of the graph, returns unevaluated.
