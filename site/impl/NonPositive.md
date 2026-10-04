---
source: src/core.c
---
**Algorithm.** `builtin_nonpositive` (`src/core.c`) decides whether a numeric
quantity is real and `<= 0`. It gates on `is_numeric_quantity` (a non-numeric
argument returns `NULL` and stays unevaluated), then classifies with the shared
`numeric_real_sign` helper: exact `EXPR_INTEGER` / `EXPR_BIGINT` / `Rational` /
`EXPR_REAL` (and `EXPR_MPFR`) read directly via `expr_numeric_sign`, everything
else `numericalize`d at machine precision. A non-real complex value answers
`False`; otherwise the verdict is `sign <= 0`, so zero is included.

**Data structures.** All `Expr`. A packed list or `NDArray` is read in one pass by
`ndint_sign_predicate(res, NDSP_NONPOSITIVE)` (`src/ndinteger.c`) into a `List` of
`True`/`False`; complex buffers, rank > 1, and `Indeterminate` elements are
declined, with `ndarray_delist_and_reeval` falling back to the scalar path that
`Listable` threading reaches for an unpacked list.

**Complexity / limits.** `O(1)` per scalar, `O(n)` over a buffer. For inexact or
symbolic-constant inputs the verdict rests on a machine-precision
numericalization. Attributes `Listable`, `Protected`.
