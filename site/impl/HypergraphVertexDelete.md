---
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hypergraph_vertex_delete` removes the named vertices **and
every hyperedge incident to one of them**, exactly as `VertexDelete` does for a
`Graph`. `edit_items` reads the second argument as one vertex or a List of
vertices; each is resolved through the memo's vertex index, and the call is left
unevaluated if a named vertex is absent. A keep-mask `kv` over vertices and a
keep-mask `ke` over hyperedges are formed (a hyperedge is dropped as soon as one
of its distinct members is deleted), and `rebuild` emits the surviving vertices
and hyperedges in their original order.

**Data structures.** The memo's vertex `GraphVIdx`, the distinct-vertex CSR
`soff/sv` (incidence test per hyperedge), and two `unsigned char` keep-masks;
`rebuild` copies the survivors into a fresh `Hypergraph`.

**Complexity / limits.** `O(n + Σ|e|)`. Mutators unshare the memoized node, so the
input object is never changed.
