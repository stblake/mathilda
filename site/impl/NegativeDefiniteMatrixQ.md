---
source: src/linalg/negdef_q.c
---
**Algorithm.** `builtin_negative_definite_matrix_q` tests `Re[Conjugate[x].m.x] < 0` for all nonzero `x`, equivalently that `−m` is positive definite. It mirrors `PositiveDefiniteMatrixQ` exactly but negates entries at load time: after the square-matrix shape gate and the `(re, im)`-double coercion (`ndq_leaf_to_double`), it loads `−m` into a column-major buffer, forms the Hermitian part of `−m` in the upper triangle, checks its diagonal is strictly positive (i.e. `m`'s diagonal strictly negative), and runs Cholesky via LAPACK `dpotrf`/`zpotrf` (with an in-house fallback). `info == 0` on `−m`'s Hermitian part ⇔ `m` negative definite ⇒ `True`; non-numeric/non-coercible entries give `False`. Wrong arity emits `NegativeDefiniteMatrixQ::argx`.
