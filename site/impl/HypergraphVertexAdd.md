---
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hypergraph_vertex_add` appends vertices to `h`, leaving
the hyperedges untouched. `edit_items` interprets the second argument: for the
Vertex heads any `List` is a list of vertices (so a List-valued vertex must be
wrapped as `{{...}}`), otherwise it is a single vertex. Each candidate not already
in the memo's vertex index — and not a duplicate among the candidates, screened by
a throwaway `GraphVIdx extra` — is appended in order after the existing vertices.
A vertex already present is silently skipped. The result is a fresh
`Hypergraph[{...}, edges]`.

**Data structures.** The memo's vertex `GraphVIdx` (membership test) and a
temporary `GraphVIdx extra` (de-duplicating the new candidates); the vertex and
(shared) edge Lists are rebuilt with `mk_hyp`/`mk_list`. Because the memo holds an
immutable reference, the original object is never mutated.

**Complexity / limits.** `O(n + k)` for `k` candidates. Only the vertex set grows;
to add hyperedges use `HypergraphEdgeAdd`.
