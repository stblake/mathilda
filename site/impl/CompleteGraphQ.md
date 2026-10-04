---
source: src/graph/graphprops.c
---
**Algorithm.** `builtin_complete_graph_q` tests whether every ordered pair of
distinct vertices `(u, v)` is joined by an edge usable from `u` to `v` — an
`UndirectedEdge` between them, or a `DirectedEdge[u, v]`, so a complete *directed*
graph needs both `u -> v` and `v -> u`. It first reads `graph_directed_edge_count`
(an `O(1)` memo query that also validates `g`); a non-graph gives that count `< 0`
and the head returns `False`. For a simple graph with a single edge kind, completeness
reduces to an edge count — `n(n−1)/2` undirected edges or `n(n−1)` directed ones, since
parallel edges cannot occur — so no traversal is needed. A mixed graph (both edge kinds)
may pair `u<->v` with `u->v`, so it falls back to the adjacency scan `induced_complete`,
which stamps the distinct out-neighbours of each vertex and checks that each reaches the
other `k−1`. Graphs on 0 or 1 vertices are complete (no pair can fail).

**Data structures.** `CompleteGraphQ[g, vlist]` builds the integer-indexed `GraphAdj`
(successor/predecessor CSR) and a `char` selection mask `sel[]` over the vertices named
in `vlist`, resolving each through `graph_vertex_position` (the memoized vertex index);
if any element of `vlist` is not a vertex of `g`, the result is `False`. `induced_complete`
uses a per-`u` `stamp[]` array so an undirected edge listed alongside a redundant directed
one is not double-counted. The graph itself is the ordinary `Graph[List verts, List edges]`
`Expr` tree; validation and the edge count come from its per-node memo.

**Complexity / limits.** The single-edge-kind fast path is `O(1)` after validation; the
induced/mixed scan is `O(V + E)`. An allocation failure in the scan leaves the call
unevaluated (`NULL`) rather than answering `False`. Every non-graph argument returns
`False`, never an unevaluated expression.
