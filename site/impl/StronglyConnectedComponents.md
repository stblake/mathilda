---
references:
  - "R. E. Tarjan, *Depth-first search and linear graph algorithms*, SIAM J. Comput. **1** (1972) 146-160."
source: src/graph/components.c
---
**Algorithm.** `builtin_strongly_connected_components` always computes strongly connected
components with an iterative Tarjan scan (`graph_strong_label`), never the undirected flood-fill
that `ConnectedComponents` falls back to. Two vertices share a component when each is reachable
from the other following edge directions; an undirected edge is traversed both ways. The
components are returned in Tarjan completion order — the same order the directed form of
`ConnectedComponents` uses.

**Data structures.** A CSR `GraphAdj` and Tarjan's `index[]`, `low[]`, `onstack[]` byte mask,
an explicit vertex stack, a per-node child cursor and an explicit recursion stack (no C
recursion). The component order array is the identity, since Tarjan's order is already the
answer order.

**Complexity / limits.** `O(V + E)`. Exactly one argument; a non-graph argument returns
unevaluated. (Mathematica spells this operation `ConnectedComponents`; Mathilda keeps the
explicit name as well.)
