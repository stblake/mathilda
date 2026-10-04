---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §§22.2 (BFS) and 24.3 (Dijkstra)."
source: src/graph/shortestpath.c
---
**Algorithm.** `GraphDistance[g, s, t]` gives the length of a shortest path from
`s` to `t`. The core `builtin_graph_distance` dispatches on
`graph_weights_usable(g)`. The unweighted default runs a breadth-first search
over the successor adjacency `GraphAdj.out[]` — for a directed graph this follows
edge direction, and for an undirected graph `out[]` is symmetric, so it is an
ordinary shortest path — and returns the integer hop count. When `g` carries a
non-negative numeric `EdgeWeight`, it runs Dijkstra over a call-scoped weighted
adjacency `WAdj` (built fresh rather than widening the shared `GraphAdj`, which
has no weight storage) and returns the accumulated distance as a machine real.
An unreachable target is `Infinity`; a symbolic or negative weight falls back to
unit BFS rather than erroring. The registered head is actually the wrapper
`builtin_gmet_graph_distance` (`gmet_distance.c`): it answers the single-source
form `GraphDistance[g, s]` (all distances from `s`) itself and forwards the
three-argument form here.

**Data structures.** The unweighted path uses `graph_build_adj`'s CSR adjacency
(`out[]`/`in[]` over one `block` allocation) with `parent[]`/`dist[]` arrays
(`-1` = unreached). The weighted path builds `WAdj` in two passes — an
`outdeg` count, then a fill — indexing endpoints through a `GraphVIdx` and owning
the resolved per-edge weight expressions; Dijkstra keeps a `done[]` flag array
and `double dist[]` (`DBL_MAX` = unreached).

**Complexity / limits.** BFS is `O(V + E)`. Dijkstra is implemented as a plain
`O(V²)` array scan (no binary heap), matching this subsystem's small-graph
exact-algorithm precedent. A weighted distance is a machine `double` (so `12.`,
not `12`), bit-for-bit consistent with `GraphDistance[g, s]` and
`GraphDistanceMatrix`. Negative edge weights are out of scope — they silently
demote the query to an unweighted hop count.
