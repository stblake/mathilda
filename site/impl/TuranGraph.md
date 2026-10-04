---
references:
  - "P. Turán, *On an extremal problem in graph theory* (Hungarian), Mat. Fiz. Lapok **48** (1941) 436-452."
source: src/graph/gmet_generators.c
---
**Algorithm.** `builtin_turan_graph` builds the Turán graph `T(n, k)`: the complete
`k`-partite graph on `n` vertices with the parts as equal as possible. The part sizes are
`sz[i] = n/k + (i < n mod k ? 1 : 0)`, so the first `n mod k` parts get one extra vertex
(larger parts first). It then emits every edge between vertices in different parts via
`multipartite`. The vertices are the integers `1..n`; `k` is clamped to `n` when larger.
This is the graph that, by Turán's theorem, has the most edges among `n`-vertex graphs with no
`(k+1)`-clique.

**Data structures.** The shared generator emits edges as packed 64-bit endpoint keys into a
`Pairs` accumulator, sorts and deduplicates them, and wraps the result as an undirected
`Graph[Range[n], {...}]`. All edges share one `UndirectedEdge` head node and the integer vertex
nodes are shared by reference.

**Complexity / limits.** `O(V + E)` to emit plus `O(E log E)` to sort when not already ordered;
`E` is the Turán edge count. Guarded by a vertex cap `GEN_MAX_VERTICES = 10^8` and an edge cap
`GEN_MAX_EDGES = 5·10^7` (either returns unevaluated). Both arguments must be positive integers;
an invalid argument returns unevaluated. Always undirected.
