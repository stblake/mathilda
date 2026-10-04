---
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_uniform_hypergraph_q` scans the raw hyperedge arities
`eoff[j+1] - eoff[j]` and returns `True` iff they are all equal. With a second
argument `k` (a non-negative integer), the common arity must equal `k`
(`k`-uniform). The first arity seen becomes the reference; a hypergraph with no
hyperedges is `True`. On a non-hypergraph the one-argument form returns `False`
(it is a `*Q` predicate); the `k` form with a malformed `k` is left unevaluated.

**Data structures.** The borrowed `HypView` from the validated-hypergraph memo;
only the raw hyperedge offsets `eoff` are read. The result is a `True`/`False`
symbol.

**Complexity / limits.** `O(m)`, short-circuiting on the first mismatch. Arity is
the `Length` as written, so a repeated vertex counts toward uniformity.
