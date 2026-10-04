---
source: src/list/distance.c
---
**Algorithm.** `builtin_cosine_distance` computes
`1 - (u . Conjugate[v]) / (Norm[u] Norm[v])`. `dist_dot_conj` forms the numerator
`Sum_i u_i Conjugate[v_i]` (the `Conjugate` makes it correct for complex vectors
and is a no-op on reals, as Mathematica writes it), and `dist_norm` gives each
Euclidean norm as `Sqrt[Sum Abs[u_i]^2]`. All of it is composed from the internal
arithmetic primitives through `eval_and_free`, so exact input stays exact and
symbolic input survives.

**Range and conventions.** The value runs over `[0, 2]`: `0` for parallel, `1`
for orthogonal, `2` for antiparallel. Unlike the Euclidean family this is *not* a
metric (it violates the triangle inequality) and has no squared form that ranks
identically, so callers use it directly. A zero vector on either side gives `0`
(a special case, since the quotient would be `0/0` → `Indeterminate`), matching
Mathematica.

**Shape / limits.** `dist_shape` admits two scalars or two equal-length `List`s;
a length mismatch or a list-valued component declines. O(n) arithmetic
evaluations, interpreter-speed. `ATTR_PROTECTED`. Shares the `dist_sum`/`dist_norm`
machinery with `EuclideanDistance`, `SquaredEuclideanDistance` and
`ManhattanDistance`.
