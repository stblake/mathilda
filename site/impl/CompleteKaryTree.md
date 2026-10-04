---
source: src/graph/gmet_generators.c
---
**Algorithm.** `builtin_complete_kary_tree` accepts `CompleteKaryTree[n]` (arity 2) or `CompleteKaryTree[n, k]`, where `n` is the number of **levels** and `k` the number of children per node. It sums the level sizes `1 + k + k^2 + ...` to get the vertex count, then hands off to `kary_tree`, which numbers vertices in heap order: vertex `i` (0-based) has children `k*i + 1 .. k*i + k`, those below the vertex count being kept.

**Data structures.** Edges go into the shared packed-64-bit `Pairs` buffer and are produced already in sorted order, so `pairs_graph` skips the sort and builds `Graph[Range[N], {UndirectedEdge[parent, child], ...}]` directly. The result is an ordinary `Expr` tree with vertices `1..N`.

**Complexity / limits.** `O(N)` for `N` vertices, with `N - 1` edges. The level-size loop checks overflow and refuses a tree exceeding 10^8 vertices; non-integer, zero or negative arguments leave the call unevaluated. `KaryTree[n, k]` shares the same builder but takes a vertex count instead of a level count.
