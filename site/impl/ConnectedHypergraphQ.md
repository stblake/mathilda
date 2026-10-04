---
source: src/graph/hyp_ops.c
---
**Algorithm.** `builtin_connected_hypergraph_q` returns `True` iff `h` has at
least one vertex and exactly one connected component. It runs the same
`vertex_uf` union–find as `HypergraphConnectedComponents` (unioning the members of
each hyperedge's distinct set) and counts the roots `p[i] == i`; the answer is
`roots == 1`. It is the Wolfram Function Repository name and accepts a bare List
of hyperedges (`hyp_arg`); a non-hypergraph argument gives `False` (a `*Q`
predicate), and the empty List `{}` gives `False` (the FR function leaves it
unevaluated).

**Data structures.** The distinct-vertex CSR `soff/sv`; a union–find parent array
(path-halving). The result is a `True`/`False` symbol.

**Complexity / limits.** Near-linear, `O(Σ|e| · α(n))` for the union step plus
`O(n)` to count roots.
