---
references:
  - "N. Chiba and T. Nishizeki, *Arboricity and subgraph listing algorithms*, SIAM J. Comput. **14** (1985) 210-223."
source: src/graph/gmet_cluster.c
---
**Algorithm.** `builtin_global_clustering_coefficient` gives
`3 x (triangles) / (connected triples)` exactly — equivalently
`sum t(v) / sum C(d(v), 2)` where `t(v)` is the number of triangles through
vertex `v` and `d(v)` its degree. It shares the triangle machinery of
`GraphTriangleCount`/`LocalClusteringCoefficient`: `tri_compute` lists triangles
over the underlying simple graph using the **degree-ordered orientation** (each
edge points from lower to higher `(degree, index)` rank), which bounds every
oriented out-degree by `O(sqrt m)` and total work by `O(m^1.5)`. For a directed
graph the denominator per vertex is `in(v) out(v) - r(v)` (reciprocal neighbours)
and a "triangle" is a directed 3-cycle. The ratio is returned in lowest terms as
an exact `Rational` (or `0` when the denominator is `0`). A mixed graph is left
unevaluated.

**Data structures.** CSR oriented adjacency (`ooff`/`oadj`, plus a per-edge
arc-flag byte `ofl` when directed), per-thread mark arrays and triangle counters,
and GMP for the final exact ratio. The listing is run on the thread team
(`gmet_parallel_for`).

**Complexity / limits.** `O(m^1.5)` time. Exact (`Integer`/`Rational`);
`EdgeWeight` is ignored. `MeanClusteringCoefficient` is the *mean of the local*
coefficients and is generally a different number from this global one.
