---
references:
  - "D. J. Watts and S. H. Strogatz, *Collective dynamics of 'small-world' networks*, Nature **393** (1998) 440-442."
  - "N. Chiba and T. Nishizeki, *Arboricity and subgraph listing algorithms*, SIAM J. Comput. **14** (1985) 210-223."
source: src/graph/gmet_cluster.c
---
**Algorithm.** `builtin_mean_clustering_coefficient` averages the local clustering coefficients. For an undirected graph the local value is `t(v) / C(d(v), 2)`, with `t(v)` the triangles through `v` and `d(v)` its degree, and `0` when `d(v) < 2`. Triangles are listed over the underlying simple graph with the degree-ordered orientation (each edge points from lower to higher `(degree, index)` rank). For a directed graph a triangle is a directed 3-cycle, and the local denominator is `in(v) out(v) - r(v)`, where `r(v)` counts neighbours joined in both directions. The directed convention was reverse-engineered from Mathematica 15. The mean is exact: numerators are grouped by denominator and summed as GMP rationals (`mpq_t`), then divided by `n`. The empty graph gives `0`.

**Data structures.** `tri_compute` fills per-vertex `int64_t` arrays `t` and `den`. Each oriented edge carries a two-bit record of which arcs exist, so the directed 3-cycle test is two bit-ANDs per triangle. Rows are processed on the thread team with per-thread mark arrays. The result is an Integer or Rational `Expr`.

**Complexity / limits.** `O(m^1.5)` total, since every oriented out-degree is bounded by `O(sqrt m)`. Weights are ignored, and a mixed graph (directed and undirected edges together) is left unevaluated.
