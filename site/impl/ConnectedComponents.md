---
references:
  - "R. E. Tarjan, *Depth-first search and linear graph algorithms*, SIAM J. Comput. **1** (1972) 146-160."
source: src/graph/components.c
---
**Algorithm.** `builtin_connected_components` dispatches on direction. If the graph has any
directed edge it computes the **strongly connected components** with an iterative Tarjan
scan (`graph_strong_label`); otherwise it labels the **undirected components** by an iterative
DFS flood-fill over the combined out+in adjacency, then stably counting-sorts the components so
the largest comes first (ties broken by first appearance). Within each component the vertices
are kept in `VertexList` order.

**Data structures.** A CSR `GraphAdj` built from the validated graph. Tarjan uses `index[]`,
`low[]`, an `onstack[]` byte mask, an explicit vertex stack, a per-node child cursor and an
explicit recursion stack (no C recursion); the undirected path uses a `comp[]` label array with
an explicit DFS stack. A `keep[]` mask supports the selection form `ConnectedComponents[g, {v,
...}]`.

**Complexity / limits.** `O(V + E)` for either labelling, plus an `O(n+k)` counting sort. No
cap. The directed (Tarjan) components come out in completion order — a reverse topological order
of the condensation, so there is no edge from a component to a later one, matching Mathematica.
A non-graph argument returns unevaluated.
