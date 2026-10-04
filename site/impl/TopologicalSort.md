---
references:
  - "A. B. Kahn, *Topological sorting of large networks*, Communications of the ACM **5** (1962) 558-562."
source: src/graph/acyclic.c
---
**Algorithm.** `builtin_topological_sort` orders the vertices of a directed
acyclic graph so that `u` precedes `v` for every edge `u -> v`, by **Kahn's
algorithm**: repeatedly remove a vertex of in-degree 0 and decrement its
successors' in-degrees. Ties — several vertices ready at once — are broken by
**VertexList position** through a binary min-heap keyed on vertex index, so the
smallest ready index always goes first; the order is therefore deterministic and
an already-sorted vertex list comes back unchanged. The call is left
**unevaluated** (returns `NULL`) for a graph that has a cycle (Kahn's queue
drains before all vertices are emitted) or that carries any **undirected** edge
(an undirected edge imposes no order, so the graph is not a DAG); an edgeless
graph sorts to its `VertexList`. `TopologicalSort[{rules}]` builds `Graph[{rules}]`
through the ordinary constructor first — so the edge sugar and validation are
exactly `Graph`'s — then sorts that; a non-graph, non-list argument declines.

**Data structures.** The adjacency comes from `graph_build_adj` as a `GraphAdj`
(forward out-neighbour lists plus per-vertex in-degrees), read straight off the
validated-graph memo. Working state is an `indeg[]` counter array and the
integer min-heap; the result is assembled as a `List` of `expr_copy`'d vertex
expressions. That list is cached on the graph node's memo slot
(`GRAPH_CACHED_TOPOSORT`): a `List` is an immutable value, so handing out
references to the one cached copy is safe, and a repeated query on the same graph
node is free.

**Complexity / limits.** Linear in the graph, `O(V + E)` plus the `O(V log V)`
heap work the deterministic tie-break costs; `O(1)` on a memo hit. The output is a
plain vertex `List`, not a graph. There is no partial-order output for a cyclic or
mixed graph — the head declines rather than return a best effort, matching
Mathematica's leaving `TopologicalSort` unevaluated there.
