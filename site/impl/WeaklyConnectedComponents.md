---
source: src/graph/components.c
---
**Algorithm.** `builtin_weakly_connected_components` labels the components of the *underlying
undirected graph* with an iterative DFS flood-fill over the combined out+in adjacency — edge
directions are ignored — then stably counting-sorts the components largest first (ties broken by
first appearance). It never runs Tarjan's strongly-connected scan. Vertices inside each
component keep `VertexList` order.

**Data structures.** A CSR `GraphAdj`, a `comp[]` label array, an explicit DFS stack, and the
size buckets of `order_by_size_desc`. A `keep[]` mask supports the selection form
`WeaklyConnectedComponents[g, {v, ...}]`.

**Complexity / limits.** `O(V + E)` for the labelling plus an `O(n+k)` counting sort. No cap. A
non-graph argument returns unevaluated.
