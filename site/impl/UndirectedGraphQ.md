---
source: src/graph/graphprops.c
---
**Algorithm.** `builtin_undirected_graph_q` returns `True` exactly when `g` has no
directed edges: it is the single test `graph_directed_edge_count(g) == 0`. An edgeless
graph is therefore undirected (there is nothing oriented), matching the Wolfram Language.
A graph that mixes directed and undirected edges, or is purely directed, gives `False`.

**Data structures.** `graph_directed_edge_count` is an `O(1)` query against `g`'s
validated-graph memo (the per-node cache described in `graph.h`): on the first call it
walks `g`'s edge list once to count `DirectedEdge` nodes and memoizes the vertex index and
edge-key set; later calls read the cached count directly. The graph is the plain
`Graph[List verts, List edges]` `Expr` tree — no separate representation is materialized.

**Complexity / limits.** `O(1)` on a memo hit, one `O(E)` pass on the first query. A
non-graph argument makes `graph_directed_edge_count` return `−1`, so the `== 0` test is
`False` — the predicate never leaves the call unevaluated.
