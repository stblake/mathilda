---
references:
  - "S. G. Aksoy, C. Joslyn, C. Ortiz Marrero, B. Praggastis and E. Purvine, *Hypernetwork science via high-order hypergraph walks*, EPJ Data Science **9**:16 (2020)."
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hypergraph_line_graph` gives the line graph on the
hyperedge indices `1..m`, with `i <-> j` when hyperedges `i` and `j` intersect; a
second argument `s` requests the **s-line graph** (`|e_i ∩ e_j| >= s`, after Aksoy
et al.). Intersections are counted through the shared `s_neighbours` helper: for
hyperedge `i`, it walks `i`'s distinct vertices and, through the vertex→hyperedge
incidence CSR, increments a per-hyperedge counter, collecting each touched `j > i`
whose overlap reaches `s`, then zeroes the scratch. Edges are emitted in `(i, j)`
order (neighbours sorted). Every hyperedge index is a vertex, even when isolated.

**Data structures.** The distinct-vertex CSR `soff/sv` and the incidence CSR
`voff/ve` (`hyp_view(..., 1)`); reusable `cnt`/`touched`/`nb` scratch arrays left
zeroed between hyperedges; an `EVec` of `UndirectedEdge` expressions. The result
is a `Graph` on the index integers `1..m`.

**Complexity / limits.** `O(Σ_v deg(v)²)` — the number of hyperedge pairs meeting
at a vertex, which bounds the 1-line graph's edge count anyway. `s` must be a
positive integer, else the call is left unevaluated.
