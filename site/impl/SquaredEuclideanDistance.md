---
source: src/list/distance.c
---
**Algorithm.** `builtin_squared_euclidean_distance` computes
`Sum_i Abs[u_i - v_i]^2` via the shared `dist_builtin(res, p=2, root=false)` —
the Euclidean distance without the final square root. The sum is assembled by
`dist_sum` from the internal arithmetic primitives through `eval_and_free`, so a
real squared term skips the redundant `Abs`, a complex term uses its modulus
(`Abs`-then-square), and a symbolic term survives.

**Why the squared form matters.** Because no root is taken, the result is
*rational for rational input* (`SquaredEuclideanDistance[{1/3, 0}, {0, 1/7}]` is
`58/441`, not a float), and squaring is monotone on non-negatives. So ranking on
the squared distance orders points identically to ranking on the true distance
without ever introducing an irrational — which is exactly what lets
`FindClusters` partition n-dimensional exact data exactly (it is the metric the
spanning-tree builder ranks on).

**Shape / limits.** `dist_shape` admits two scalars or two equal-length `List`s;
a length mismatch or a list-valued component declines. O(n) arithmetic
evaluations, interpreter-speed. `ATTR_PROTECTED`. See `EuclideanDistance` for the
rooted form and `ManhattanDistance` / `CosineDistance` for the siblings.
