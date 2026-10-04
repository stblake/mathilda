---
source: src/graph/gops_setops.c
---
**Algorithm.** `builtin_graph_difference` takes exactly two graphs and builds the graph on the
union of their vertex sets whose edges are those of `g1` that are *not* in `g2`. Vertices are
merged in canonical (`expr_compare`) order. Each edge of `g1` is tested for membership in `g2`
(mapping endpoints onto the shared vertex index) and kept when absent; the kept edges are
counting-sorted into canonical order. The orientation of `g1`'s edges is preserved and edge
weights are dropped.

**Data structures.** A `GopsView` per graph with integer endpoints and a `GopsKeySet` hash of
`g1`'s edges probed by `g2`'s; the vertex-union sort is radix/counting for all-integer vertices
and a comparison sort otherwise.

**Complexity / limits.** `O(V + E)` apart from the `O(V log V)` vertex-union sort. Exactly two
arguments; a non-graph argument returns unevaluated. An undirected edge equals its reversal under
the edge key, so orientation does not affect the difference.
