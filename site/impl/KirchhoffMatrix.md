---
source: src/graph/gmet_distance.c
---
**Algorithm.** `builtin_kirchhoff_matrix` forms the Laplacian `D - A` in a single pass over the edge list. For each edge `k` with endpoints `a = eu[k]`, `b = ev[k]` it increments both diagonal entries, so each vertex's diagonal entry is the number of edges incident to it. It then sets `K[a][b] = -1`, and also `K[b][a] = -1` when the edge is undirected. Edge weights are ignored. A directed edge `u -> v` therefore contributes to both degrees but only to the `u`-to-`v` off-diagonal entry.

**Data structures.** The edge-index view from `graph_edge_indices` supplies the parallel arrays `eu`, `ev` and `directed`. The result is a dense `n*n` `int64_t` buffer returned as a packed int64 `NDArray` (`ndbuild_open`), with a plain nested `List` of integers as the fallback. Mathematica returns a `SparseArray`; Mathilda has none, so this is its `Normal` form.

**Complexity / limits.** `O(n^2)` time and space for the zero-fill, plus `O(m)` for the scan. Repeated edges accumulate on the diagonal but not off it. The empty graph and non-graph arguments leave the call unevaluated.
