---
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hypergraph_rank` is `rank_impl(res, want_max = 1)`: over
the raw hyperedge CSR `eoff` it tracks the largest arity `eoff[j+1] - eoff[j]`,
counting a repeated vertex (the `Length`, not the distinct-vertex count). With no
hyperedges it is `0`. `HypergraphCorank` is the same routine with
`want_max = 0`.

**Data structures.** The borrowed `HypView` from the validated-hypergraph memo;
the raw hyperedge offsets `eoff` are all that is read. The result is a single
`Integer`.

**Complexity / limits.** `O(m)` — one pass over the hyperedge offsets. Left
unevaluated on a non-hypergraph.
