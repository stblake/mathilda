---
references:
  - "S. G. Aksoy, C. Joslyn, C. Ortiz Marrero, B. Praggastis and E. Purvine, *Hypernetwork science via high-order hypergraph walks*, EPJ Data Science **9**:16 (2020)."
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hyperedge_connected_components` groups the hyperedges
under intersection, returning Lists of 1-based hyperedge indices. The default
`s = 1` runs union–find directly over the vertex→hyperedge incidence lists
(hyperedges sharing a vertex are unioned). For `s >= 2` (the s-connectivity of
Aksoy et al.) each hyperedge is joined to the `s_neighbours` it shares at least
`s` vertices with, then the classes are closed transitively by union–find. A
hyperedge with fewer than `s` distinct vertices lies on no s-walk and is marked
inactive, so it is omitted (empty hyperedges never appear). `groups_to_list`
orders classes by smallest index.

**Data structures.** The distinct-vertex CSR `soff/sv` and incidence CSR
`voff/ve`; a union–find parent array; an `active` mask gating hyperedges by arity;
for `s >= 2` the `cnt`/`touched`/`nb` overlap scratch. The result is a `List` of
index Lists.

**Complexity / limits.** `s = 1` is near-linear in the total incidence; `s >= 2`
is `O(Σ_v deg(v)²)` (overlap counting). `s` must be a positive integer, else the
call is left unevaluated.
