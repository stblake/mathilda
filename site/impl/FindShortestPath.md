---
references:
  - "E. W. Dijkstra, *A note on two problems in connexion with graphs*, Numer. Math. **1** (1959) 269-271."
source: src/graph/shortestpath.c
---
**Algorithm.** `builtin_find_shortest_path` returns a shortest `s`-`t` path as a list of
vertices. It dispatches on `graph_weights_usable(g)`: a graph carrying an `EdgeWeight` list in
which every weight is a non-negative, non-complex number runs **Dijkstra**; otherwise — an
unweighted graph, or one with any symbolic, negative or complex weight — it falls back to
**BFS** rather than erroring. Both follow edge direction on a directed graph and treat an
undirected edge as usable both ways, and both reconstruct the path by walking a `parent[]` array
back from `t`.

**Data structures.** BFS uses the shared CSR `GraphAdj` with an integer FIFO queue, a `dist[]`
array (`-1` = unreached) and `parent[]`. Dijkstra builds a *separate* call-scoped weighted
adjacency `WAdj` (via a `GraphVIdx` hash) rather than widening the shared `GraphAdj`, with a
`double` `dist[]` (`DBL_MAX` = unreached), a `done[]` flag array and `parent[]`; it selects the
next vertex by a linear array scan, not a heap.

**Complexity / limits.** BFS is `O(V + E)`; the array-scan Dijkstra is `O(V^2)`. When there is
no `s`-`t` path the result is `{}` (the empty list). A non-graph argument, or an `s`/`t` that is
not a vertex, returns unevaluated. (The companion `GraphDistance` returns the path *length*, or
`Infinity` when unreachable.)
