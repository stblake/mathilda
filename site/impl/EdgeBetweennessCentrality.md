---
references:
  - "U. Brandes, *A faster algorithm for betweenness centrality*, J. Math. Sociology **25** (2001) 163-177."
source: src/graph/gmet_centrality.c
---
**Algorithm.** `builtin_edge_betweenness_centrality` is the edge variant of **Brandes'
algorithm**: the same reverse-order dependency accumulation, but each unit of dependency is
charged to `acc[eid]` — the id of the arc carrying it — instead of to a vertex. Unlike the
vertex form, it *does* use `EdgeWeight` as edge lengths, so when the graph carries usable
weights it takes the Dijkstra variant (`gmet_csr_build(..., ew, want_eid=1, ...)`), and with
equal path lengths within a relative tolerance `BRANDES_TIE = 1e-12` it splits the flow across
the tied shortest paths. It sums over ordered pairs on undirected graphs too (no ×0.5),
matching the convention of the reference implementation.

**Data structures.** The per-source `GmetCSR` carries an `eid[]` array mapping each arc back to
its `EdgeList` position; the accumulator has one slot per edge rather than per vertex. The rest
is shared with the vertex form: `dist[]`/`sigma[]`/`delta[]` and the settled-vertex stack, run
per source across a thread team and merged. The result is a packed machine-real vector in
`EdgeList` order.

**Complexity / limits.** `O(nm)` unweighted, `O(nm + n^2 log n)` weighted, `n = |V|`,
`m = |E|`. No size cap; non-graph input returns unevaluated.
