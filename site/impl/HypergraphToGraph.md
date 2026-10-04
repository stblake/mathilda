---
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hypergraph_to_graph` reads `h` as an **ordered**
hypergraph (the Wolfram Function Repository convention): hyperedge
`{v1, ..., vk}` contributes the directed edge `v_a -> v_b` for every position
`a < b`. It therefore reads the **raw** hyperedge CSR `ev` (order and repeats
preserved), not the distinct set. Because Mathilda graphs are simple, a self-loop
(from a repeated vertex, `x == y`) is dropped and parallel copies are merged via a
`U64Set` on `pair_key(x, y)`; every vertex of `h` is kept, including one that
occurs only in a unary hyperedge (the FR function drops those). Accepts a
`Hypergraph` or a bare List of hyperedges (`hyp_arg`).

**Data structures.** The raw hyperedge CSR `eoff/ev`; a `U64Set` for
directed-edge de-duplication; an `EVec` of `DirectedEdge` expressions. The result
is a `Graph` on `VertexList[h]`.

**Complexity / limits.** `O(Σ_j |e_j|²)` directed pairs, de-duplicated in
amortised `O(1)` each. The simple-graph normalisation (dropped self-loops, merged
parallels, retained vertices) is the documented deviation from the FR multigraph.
