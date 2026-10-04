---
source: src/graph/hyp_util.c
---
**Algorithm.** `builtin_hypergraph_q` is a one-argument predicate that returns
`True`/`False` from `hypergraph_is_valid(h)`. That first checks the head is the
interned `Hypergraph` symbol, then calls `hyp_memo(h)`, which returns a memo slot
(validating and memoizing on a miss) or `NULL`. Validation requires the
`Hypergraph[List, List]` shape, pairwise-distinct vertices (a repeated vertex in
the vertex List is not canonical), and every hyperedge a `List` of declared
vertices. A `Graph`, a bare List of hyperedges, and a malformed — hence
unevaluated — `Hypergraph[...]` all give `False`; correspondingly `GraphQ` of a
hypergraph is `False`.

**Data structures.** No result tree beyond the `True`/`False` symbol; the work is
the shared validated-hypergraph memo (vertex `GraphVIdx` index plus the
vertex-index CSR of the hyperedges) described on the `Hypergraph` page.

**Complexity / limits.** `O(1)` on a memoized object (pointer-keyed slot hit);
`O(Σ|e|)` on the first validation of a new object, which also fills the memo. The
predicate has no side effects and never mutates its argument.
