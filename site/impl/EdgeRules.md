---
source: src/graph/gops_edit.c
---
**Algorithm.** `builtin_edge_rules` maps each edge of the graph to a `Rule` `u -> v`, in
`EdgeList` order. It reads the edge list directly and, for each edge, copies its two endpoints
and wraps them in `Rule` — regardless of whether the edge was a `DirectedEdge` or an
`UndirectedEdge`. Direction is therefore *not* encoded in the output: every edge becomes a plain
`->` rule, which is the form the `Graph` constructor accepts as input sugar.

**Data structures.** None beyond the result `List` of `Rule` nodes; no adjacency or vertex index
is built, and the edge order is the graph's stored order.

**Complexity / limits.** `O(E)`. Exactly one argument; a non-graph argument returns unevaluated.
The result is a plain list of rules, not a `Graph`.
