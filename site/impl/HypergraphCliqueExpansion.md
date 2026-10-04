---
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_hypergraph_clique_expansion` gives the 2-section: the
undirected `Graph` on `VertexList[h]` with `u <-> v` whenever `u` and `v` share a
hyperedge. Each hyperedge is read as its distinct-vertex set `sv[soff[j]..]`; for
every pair `(a, b)` with `a < b` in that set, the oriented key `pair_key(min, max)`
is inserted into an open-addressing `U64Set`, and a first insertion emits an
`UndirectedEdge` in first-co-occurrence order. A vertex repeated inside a
hyperedge yields no self-loop (the distinct set collapses it). The argument may be
a `Hypergraph` or, for Function-Repository compatibility, a bare List of
hyperedges (`hyp_arg` wraps and evaluates it).

**Data structures.** The distinct-vertex CSR `soff/sv` from the memo; a growable
`U64Set` (64-bit packed-pair keys, multiplicative hash) for edge de-duplication;
a growable `EVec` of `UndirectedEdge` expressions. The result is a `Graph`, itself
validated and memoized by the Graph subsystem.

**Complexity / limits.** `O(Σ_v deg(v)²)` in the worst case — the number of
co-occurring pairs, which bounds the output edge count anyway. This is also the
graph whose shortest-path distance `HypergraphDistance` measures.
