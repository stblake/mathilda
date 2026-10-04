---
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hypergraph_connected_components` gives the vertex
components — two vertices are connected when a chain of hyperedges joins them.
`vertex_uf` runs union–find over the vertices: within each hyperedge's distinct
set it unions the first member with every other. `groups_to_list` then collects
the classes, ordered by their smallest member (first vertex in `VertexList`
order), vertices within a class ascending. An isolated vertex is a singleton
component. Accepts a `Hypergraph` or a bare List of hyperedges (`hyp_arg`).

**Data structures.** The distinct-vertex CSR `soff/sv`; a union–find parent array
with path-halving and union-by-smaller-root; `groups_to_list` scratch (`id`,
`size`, `cid`, per-class member buffers). The result is a `List` of vertex Lists.

**Complexity / limits.** Near-linear, `O(Σ|e| · α(n))`. Note the ordering differs
from Graph `ConnectedComponents` (which is largest-first / strong components); the
hypergraph form is deterministically ordered by first vertex.
