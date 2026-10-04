---
references:
  - "S. G. Aksoy, C. Joslyn, C. Ortiz Marrero, B. Praggastis and E. Purvine, *Hypernetwork science via high-order hypergraph walks*, EPJ Data Science **9**:16 (2020)."
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hyperedge_distance` gives the length of a shortest walk of
intersecting hyperedges from `i` to `j` (an s-walk, consecutive hyperedges sharing
at least `s` vertices, when a fourth argument `s` is given); `Infinity` if none,
`0` for `i == j`. `edge_bfs` runs a breadth-first search directly on the incidence
structure — the line graph is never materialised. For `s = 1` it expands each
vertex's incidence list at most once per BFS (a `vexp` stamp), which keeps it
linear; for `s >= 2` it expands neighbours through `s_neighbours`. The one-source
forms `HyperedgeDistance[h, i]` and `[h, i, All, s]` return the whole distance
vector (stopping early only when a specific target is named).

**Data structures.** The distinct-vertex CSR `soff/sv` and incidence CSR
`voff/ve`; a BFS queue and a `vexp` vertex-expanded bitmap (`s = 1`) or
`cnt`/`touched`/`nb` overlap scratch (`s >= 2`); `dist_list` packs the all-targets
answer into an `int64` buffer when every distance is finite, else emits a `List`
with `Infinity` entries.

**Complexity / limits.** `s = 1` is `O(Σ|e|)`; `s >= 2` is `O(Σ_v deg(v)²)`. An
out-of-range hyperedge index leaves the call unevaluated; `s` must be a positive
integer.
