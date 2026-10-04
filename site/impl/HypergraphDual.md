---
references:
  - "C. Berge, *Hypergraphs: Combinatorics of Finite Sets*, North-Holland Mathematical Library 45 (Elsevier, 1989), ch. 1 (dual hypergraph)."
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hypergraph_dual` builds the dual: vertices `1..m`, one per
hyperedge of `h`, and one hyperedge per vertex of `h`. It requests the incidence
CSR (`hyp_view(..., want_incidence = 1)`) so that, for vertex `i` in `VertexList`
order, the dual hyperedge is the list of hyperedge indices `ve[voff[i]..voff[i+1])`
— already ascending and de-duplicated, since the incidence is built from the
distinct-vertex sets. An isolated vertex gives an empty hyperedge. The result is
`Hypergraph[{1, ..., m}, {...}]`, so the dual of the dual recovers the incidence
structure on `1..n`.

**Data structures.** The incidence CSR `voff/ve` from the memo. The `1..m` vertex
`Integer`s are built once and shared (via `expr_copy`) into both the dual's vertex
List and its hyperedges; the dual is assembled with `mk_hyp`/`mk_list`.

**Complexity / limits.** `O(Σ|e|)` — linear in the total incidence, which equals
the total size of the dual's hyperedges. A bare List of hyperedges is not accepted
(the head needs the memoized incidence), so it is left unevaluated.
