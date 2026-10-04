---
references:
  - "F. Harary, *The maximum connectivity of a graph*, Proc. Natl. Acad. Sci. USA **48** (1962) 1142-1146."
source: src/graph/gmet_generators.c
---
**Algorithm.** `builtin_harary_graph` constructs `HararyGraph[k, n]`, the
`k`-connected graph on `n` vertices with the fewest possible edges (Harary 1962).
With `r = ⌊k/2⌋` the core is the circulant `C_n(1, 2, …, r)`: each vertex joined to
its `r` nearest neighbours on each side of the cycle `0, 1, …, n−1`
(`circulant_add`). When `k` is odd two cases add the last ring of edges:
`n` even adds the `n/2` "diameter" chords `i — i+n/2`; `n` odd adds a near-perfect
matching across the half-cycle (`0 — h`, `0 — h+1`, and `i — i+h+1` for
`i = 1 … h−1`, with `h = (n−1)/2`). Vertices are emitted to `pairs_graph`, which
renumbers them to `1 … n` and deduplicates parallel edges, yielding an undirected
graph that is exactly `k`-vertex-connected with `⌈kn/2⌉` edges.

**Data structures.** A `Pairs` edge accumulator (a growable array of endpoint
pairs) filled by `circulant_add`/`pairs_add`, then turned into the `Graph`
expression by `pairs_graph`, which sorts and dedups so the doubled edges of small
cases collapse.

**Complexity / limits.** `O(kn)` edges are emitted and the build is linear in that
size. The arguments must satisfy `k ≥ 0` and `n > k` (a `k`-connected graph needs
at least `k+1` vertices); `n` is capped at `GEN_MAX_VERTICES`. Anything outside
those bounds leaves the call unevaluated. Note the argument order is
`HararyGraph[k, n]` — connectivity first, vertex count second.
