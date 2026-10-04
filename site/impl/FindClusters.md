---
source: src/list/find_clusters.c
---
**Algorithm.** `builtin_find_clusters` partitions a list into clusters of nearby
elements. `fc_probe_shape` decides the element kind: real scalars (distance
`|a - b|` on the line), equal-length numeric vectors (squared Euclidean
distance), colours whose arguments are coordinates (`RGBColor`/`GrayLevel`/
`Hue`/`CMYKColor`, one head throughout), or strings (`EditDistance`); a mixture,
a ragged shape, a non-real component, `Complex`, or a visible `NDArray` declines.
The count is `Automatic` (cut a sorted-adjacent gap exceeding `FC_GAP_FACTOR = 3`
times the median gap), `UpTo[n]` (bounded), or `n` (fixed, capped at the distinct
count). Clusters appear by first occurrence; elements keep input order within a
cluster.

**Methods and exactness.** Named methods include `Agglomerate`/`SpanningTree`,
`KMeans`, `KMedoids`, `DBSCAN`, `MeanShift`, `NeighborhoodContraction`,
`JarvisPatrick`, `GaussianMixture` (the `src/ml` EM fit) and `Spectral`. In one
dimension the spanning tree *is* the sorted adjacency chain, computed on the
elements themselves via `list_numeric_cmp`, so exact 1-D input is ordered
exactly; above one dimension it is a real minimum spanning tree built by Prim's
algorithm over exact distances, with a machine-double Prim for points whose every
coordinate is already machine. The inherently-inexact methods (KMeans, the
density family, GaussianMixture, Spectral) work on a double projection, which is
correct — a mean, kernel or eigenvector is inexact by definition.

**Complexity / limits.** 1-D is O(n log n) (dominated by the sort) with no size
cap (~2.3 s at 10⁶). Multi-dimensional MST is quadratic: machine points capped at
20000, exact points and strings at 2000; `Spectral` builds an n×n matrix and
declines above `FC_SPECTRAL_MAX_N = 2000`. The result is **not** bit-identical to
Mathematica (which auto-selects an unpublished metric/preprocessing); the
`tests/test_list.c` acceptance table is the specification. `ATTR_PROTECTED`.
