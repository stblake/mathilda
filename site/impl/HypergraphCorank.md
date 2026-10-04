---
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hypergraph_corank` is `rank_impl(res, want_max = 0)`, the
mirror of `HypergraphRank`: over the raw hyperedge CSR `eoff` it keeps the
smallest arity `eoff[j+1] - eoff[j]` (the `Length`, counting a repeated vertex).
With no hyperedges it is `0`.

**Data structures.** The borrowed `HypView` from the validated-hypergraph memo;
only the raw hyperedge offsets `eoff` are read. The result is a single `Integer`.

**Complexity / limits.** `O(m)` — one pass over the hyperedge offsets. Left
unevaluated on a non-hypergraph.
