---
references:
  - "C. Berge, *Hypergraphs: Combinatorics of Finite Sets*, North-Holland Mathematical Library 45 (Elsevier, 1989), ch. 1 (induced sub-hypergraph)."
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hypergraph_restriction` gives Berge's induced
sub-hypergraph: every hyperedge of `h` intersected with the given vertex set, with
hyperedges that miss it entirely dropped. `vertex_mask` builds the membership
bitmap `in`. Each hyperedge is walked over its **raw** elements `ev` (so order and
repeats of the kept vertices survive); if all its elements are kept the original
`List` is shared unchanged, otherwise a new `List` of the kept elements is built.
An intersection that is empty is dropped. The kept vertices (in `VertexList` order)
form the new vertex set.

**Data structures.** The raw hyperedge CSR `eoff/ev` and distinct set are read
through the `HypView`; an `in` bitmap, and `vs`/`es`/`tmp` buffers for the rebuilt
vertices, hyperedges, and per-hyperedge kept elements. The result is a fresh
`Hypergraph`.

**Complexity / limits.** `O(n + Σ|e|)`. Unlike `Subhypergraph`, a hyperedge
straddling the boundary is kept in truncated form, so the result may have more
hyperedges than `Subhypergraph` on the same vertex set. A bare List is not
accepted.
