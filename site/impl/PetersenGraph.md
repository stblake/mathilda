---
references:
  - "M. E. Watkins, *A theorem on Tait colorings with an application to the generalized Petersen graphs*, J. Combin. Theory **6** (1969) 152-164."
source: src/graph/gmet_generators.c
---
**Algorithm.** `builtin_petersen_graph` emits the generalized Petersen graph
`GP(n, k)` on `2n` vertices. `PetersenGraph[]` is the default `GP(5, 2)` — the
Petersen graph itself. For `PetersenGraph[n, k]` the three edge families are added
with `pairs_add` over `i = 0 … n−1`: the inner "star" `i — (i+k) mod n`, the outer
cycle `(n+i) — (n + (i+1) mod n)`, and the spokes `i — (n+i)`. The result is handed
to `pairs_graph`, which renumbers the vertices to `1 … 2n` (inner `1 … n`, outer
`n+1 … 2n`) and deduplicates, so the doubled inner/outer edges at `n = 2` collapse
to a simple graph.

**Data structures.** A `Pairs` endpoint-pair accumulator filled in one pass, then
converted to the `Graph` expression by `pairs_graph` (which sorts and dedups the
edge list).

**Complexity / limits.** `O(n)` edges emitted and a linear build. The arguments
must satisfy `n ≥ 0`, `k ≥ 0`, `k` not a multiple of `n` (`k mod n ≠ 0`, so the
inner edges form a genuine circulant and not self-loops), and `n ≤
GEN_MAX_VERTICES/2`; otherwise the call is left unevaluated. The standard Petersen
graph `GP(5, 2)` is 3-regular on 10 vertices with 15 edges.
