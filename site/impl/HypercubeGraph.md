---
source: src/graph/gmet_generators.c
---
**Algorithm.** `builtin_hypercube_graph[d]` builds the `d`-dimensional cube graph on `2^d` vertices. Vertex numbers are `1 + v` for bit patterns `v`; for each `v` and each bit `b` that is clear in `v`, the edge `v ~ v | (1 << b)` is emitted, so two vertices are adjacent exactly when their 0-based labels differ in one bit.

**Data structures.** Edges are written to the packed-64-bit `Pairs` buffer (one `uint64_t` key per edge) and `pairs_graph` sorts them lexicographically and constructs `Graph[Range[2^d], {UndirectedEdge[i, j], ...}]` as an ordinary `Expr` tree.

**Complexity / limits.** `O(d * 2^d)` time and `d * 2^(d-1)` edges. The dimension must be a machine integer in `0..24` (`HypercubeGraph[0]` is a single vertex); larger or non-integer arguments, and anything exceeding the 5 x 10^7 edge cap, are left unevaluated.
