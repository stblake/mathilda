---
source: src/linalg/cross.c
---
**Algorithm.** `builtin_cross` implements the generalised cross product. Given `m = n-1` input vectors of length `n` (each must be a `List` of equal length, exactly one less than the number of arguments — otherwise it emits `Cross::nonn1` and returns `NULL`), the `i`-th component of the result is the signed minor obtained by stacking the input rows and deleting column `i`. Each minor is evaluated by `laplace_det` (shared with `Det`), and a sign `(-1)^(m+i)` is applied by wrapping the determinant in `Times[-1, …]` and calling `eval_and_free`. The classic 3-vector case `Cross[u, v]` is the `m = 2`, `n = 3` instance.

**Data structures.** Each minor is assembled into a flat `Expr**` array of `m*m` element pointers (borrowed from the input vectors, not copied) and passed to `laplace_det` with an explicit column-index list. Entries flow through symbolically, so the result is exact/symbolic whenever the inputs are. The output is a `List` of `n` components.
