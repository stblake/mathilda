---
source: src/graph/gmet_generators.c
---
**Algorithm.** `builtin_grid_graph` accepts one argument: an integer `n` (a path of `n` vertices) or a list `{n1, ..., nk}` of up to 32 positive dimensions. Vertices are numbered `1 + x1 + n1*x2 + n1*n2*x3 + ...`, so the first coordinate varies fastest. For each vertex and each axis `i` with room to step (`x_i + 1 < n_i`), it emits the edge to the vertex one stride further along that axis.

**Data structures.** Per-axis strides are held in small `int64_t[32]` arrays and edges go into the packed-64-bit `Pairs` buffer, which is already in sorted order for this loop. `pairs_graph` then returns the canonical `Graph[Range[N], {UndirectedEdge[u, v], ...}]` tree.

**Complexity / limits.** `O(N * k)` for `N = n1...nk` vertices; the exact edge count `sum_i (n_i - 1) N / n_i` is computed first so an over-cap grid is refused before any allocation. More than 10^8 vertices, 5 x 10^7 edges or 32 dimensions, and any non-positive or non-integer dimension, leave the call unevaluated.
