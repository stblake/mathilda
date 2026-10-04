---
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_subhypergraph` follows `Subgraph`'s rule: keep the named
vertices that are in `h` (in `VertexList` order; absent names ignored) and the
hyperedges lying **entirely** among them. `vertex_mask` turns the vertex-List
argument into a membership bitmap `in` over `h`'s vertices (non-vertices ignored).
A hyperedge is kept only if every vertex of its distinct set is in the mask; then
`rebuild` emits the kept vertices and hyperedges in original order.

**Data structures.** The distinct-vertex CSR `soff/sv`, the memo's vertex
`GraphVIdx` (for `vertex_mask`), an `in` bitmap and an `ke` hyperedge keep-mask;
`rebuild` constructs the fresh `Hypergraph`.

**Complexity / limits.** `O(n + Σ|e|)`. Contrast `HypergraphRestriction`, which
keeps *partial* hyperedges by intersecting rather than requiring full containment.
A bare List of hyperedges is not accepted (left unevaluated).
