---
references:
  - "N. Chiba and T. Nishizeki, *Arboricity and subgraph listing algorithms*, SIAM J. Comput. **14** (1985) 210-223."
source: src/graph/gmet_cluster.c
---
**Algorithm.** `builtin_local_clustering_coefficient` gives, per vertex, the
fraction of pairs of its neighbours that are adjacent: `t(v) / C(d(v), 2)`
(`0` when `d(v) < 2`), where `t(v)` is the number of triangles through `v`. For a
directed graph it counts directed 3-cycles through `v` over the `in(v) out(v) -
r(v)` in/out neighbour pairs. `LocalClusteringCoefficient[g]` returns the whole
list in `VertexList` order; `[g, v]` the single value. The per-vertex triangle
counts come from the same degree-ordered triangle listing as
`GraphTriangleCount`, run with `per_vertex = 1`.

**Data structures.** The CSR oriented adjacency and thread-team triangle listing
of `tri_compute` (`TriInfo.t[v]`, `TriInfo.den[v]`). Because degree sequences
repeat, most coefficients repeat, so the whole-graph form builds each distinct
`(t, den)` value once and shares it by reference through a small open-addressed
hash; the result is cached on the graph node.

**Complexity / limits.** `O(m^1.5)` for the listing, `O(V)` to assemble. Values
are exact `Rational`s (or `0`). `EdgeWeight` is ignored; a mixed graph is left
unevaluated.
