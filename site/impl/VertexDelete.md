---
source: src/graph/gops_edit.c
---
**Algorithm.** `builtin_vertex_delete` removes one or more vertices together with every incident
edge. It starts with all vertices kept, then clears the keep-flag of each deleted vertex: a
single literal vertex, each item of a list (a literal vertex, else a pattern matched against all
vertices), or a bare pattern matched against all vertices. An edge survives only when *both* of
its endpoints survive. The surviving vertices and edges are copied into a fresh canonical
`Graph`, preserving the original order and remapping endpoint indices.

**Data structures.** A `GopsView` over the raw arrays with integer endpoints, a `vkeep[]` vertex
mask and a derived `ekeep[]` edge mask (`ekeep[k] = vkeep[eu[k]] && vkeep[ev[k]]`). The view is
taken after pattern matching, which can evaluate arbitrary user code.

**Complexity / limits.** `O(V + E)`; the pattern path adds `O(V)` matches per pattern. The result
is a canonical `Graph`. A list item that is neither a vertex of the graph nor a pattern, or a
non-graph argument, returns unevaluated.
