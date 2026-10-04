---
source: src/graph/gmet_distance.c
---
**Algorithm.** `builtin_graph_distance_matrix` accepts `[g]` or `[g, d]`, where `d` is a non-negative cutoff and entries beyond it become `Infinity`. For an unweighted graph it builds a CSR of the reversed arcs and runs a bit-parallel multi-source BFS (`gmet_msbfs_run`): 256 targets per adjacency sweep, each visited vertex writing its level into a contiguous run of its row. A graph with `EdgeWeight` instead runs a binary-heap Dijkstra from every source and returns machine reals. Distances follow edge direction; an undirected edge is usable both ways. The diagonal is `0`, unreached pairs are `Infinity`.

**Data structures.** The matrix is built in a flat `n*n` `int64_t` (unweighted, `-1` for unreached) or `double` buffer. When every entry is finite it is returned directly as a packed int64 or float64 `NDArray`, with no boxing. Otherwise `matrix_from_int` / `matrix_from_real` emit a `List` of rows, each offered to the packer, with the `Infinity` symbol in the unreached slots. Source batches and per-thread BFS or heap workspaces are spread over the thread team, and results are cached by expression.

**Complexity / limits.** `O(n (n + m) / 256)` word operations unweighted and `O(n (m + n) log n)` weighted, with `O(n^2)` output. The empty graph and symbolic, complex or negative weights leave the call unevaluated.
