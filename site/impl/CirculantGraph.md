---
source: src/graph/gmet_generators.c
---
**Algorithm.** `builtin_circulant_graph` takes exactly two arguments: a positive machine-integer vertex count `n` and either one integer offset `j` or a `List` of integer offsets. Each offset is reduced into `0..n-1` (`j mod n`, negatives wrapped), and for every vertex `i` and every non-zero offset `s` the edge `i ~ (i + s) mod n` is emitted. An offset of `0` (or a multiple of `n`) contributes nothing, and an offset `s` and its mirror `n - s` give the same edges.

**Data structures.** Edges are accumulated as packed 64-bit keys `(min << 32) | max` in the shared `Pairs` buffer (`pairs_add`). `pairs_graph` sorts them lexicographically, drops duplicates (for example `s = n/2` met from both ends), and builds the canonical `Graph[Range[n], {UndirectedEdge[i, j], ...}]` expression tree. All edges share one `UndirectedEdge` head node and reuse the vertex integer nodes by reference.

**Complexity / limits.** `O(n * |offsets|)` edge generation plus an `O(E log E)` sort. Calls with a non-integer `n` or offset, a non-positive `n`, `n` above 10^8 vertices, or more than 5 x 10^7 edges are left unevaluated. Options are not supported; the result is always undirected.
