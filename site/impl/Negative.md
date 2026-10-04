---
source: src/core.c
---
**Algorithm.** `builtin_negative` (`src/core.c`) is the sign-mirror of
`builtin_positive`. It gates on `is_numeric_quantity` — a non-numeric argument
returns `NULL` and stays unevaluated — then calls the shared `numeric_real_sign`
helper, which reads exact `EXPR_INTEGER` / `EXPR_BIGINT` / `Rational` /
`EXPR_REAL` (and `EXPR_MPFR`) directly via `expr_numeric_sign` and
`numericalize`s anything else at machine precision. A value that numericalizes to
a non-real complex number answers `False`; otherwise the verdict is `sign < 0`.

**Data structures.** All `Expr`. A packed list or `NDArray` is read in one pass by
`ndint_sign_predicate(res, NDSP_NEGATIVE)` (`src/ndinteger.c`), emitting a `List`
of `True`/`False` straight off the buffer; it declines complex buffers, rank > 1,
and `Indeterminate` elements, after which `ndarray_delist_and_reeval` falls back
to the scalar path that `Listable` threading reaches for an unpacked list.

**Complexity / limits.** `O(1)` per scalar, `O(n)` over a buffer. For inexact or
symbolic-constant inputs the decision rests on a machine-precision
numericalization, so a value indistinguishable from zero at machine precision is
classified by its machine sign. Attributes `Listable`, `Protected`.
