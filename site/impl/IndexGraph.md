---
source: src/graph/gops_edit.c
---
**Algorithm.** `builtin_index_graph` handles `IndexGraph[g]` and `IndexGraph[g, r]` with `r` a machine integer (default 1). Vertex `i` of `VertexList[g]` is renamed to the integer `r + i`, so the vertices become `r, r+1, ...` in their existing order. Every edge is rebuilt from the same endpoint indices with its own directedness, so orientation and edge order are preserved, and any `EdgeWeight` list is copied across. An offset that would overflow `int64` leaves the call unevaluated.

**Data structures.** `Graph[List, List]` expression tree; because the renaming is positional, the endpoint arrays `eu`/`ev`/`directed` of the memo are reused unchanged, and only the vertex and edge expressions are new.

**Complexity / limits.** `O(V + E)`, with no hashing needed. A non-integer `r` or an invalid graph is left unevaluated.
