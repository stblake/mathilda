---
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hypergraph_edge_delete` removes every hyperedge that is
`SameQ` to a named one — compared **as written**, so `{2, 1}` does not delete
`{1, 2}`, and all copies of a repeated hyperedge are removed. `edit_items`
(`edges = 1`) reads the argument as one hyperedge or a list of hyperedges. The
named hyperedges are put into a `GraphVIdx del` keyed on the whole `List` node;
each of `h`'s hyperedges is looked up there, its keep-bit set to miss, and a `hit`
bit recorded. If some named hyperedge never matched (its first slot's `hit` is
unset) the call is left unevaluated; otherwise `rebuild` emits the survivors.

**Data structures.** A `GraphVIdx del` over the named hyperedge Lists, a `hit`
mask, and vertex/hyperedge keep-masks (vertices are all kept); `rebuild` copies the
surviving hyperedges into a fresh `Hypergraph`.

**Complexity / limits.** `O(m + k)` hyperedge comparisons (each an `expr_hash` /
`SameQ` on the List). The vertex set is unchanged — deleting a hyperedge does not
remove vertices that became isolated.
