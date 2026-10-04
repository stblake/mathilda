---
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hyperedge_sizes` takes a `HypView` of the hypergraph
(incidence not needed) and reads the **raw** hyperedge CSR `eoff`: the arity of
hyperedge `j` is `eoff[j+1] - eoff[j]`, its `Length` exactly as written, so a
repeated vertex counts and an empty hyperedge is `0`. The arities are returned in
`EdgeList` order.

**Data structures.** The borrowed integer view from the validated-hypergraph memo
(no hashing, no copy of the tree). The result is assembled through `hyp_int_list`,
which hands back a packed `int64` buffer above the packing threshold and a plain
`List` of integers otherwise.

**Complexity / limits.** `O(m)` for `m` hyperedges — the offset differences are a
single pass, independent of the total incidence. Left unevaluated on a
non-hypergraph.
