---
source: src/core.c
---
**Algorithm.** `builtin_positive` (`src/core.c`) decides the sign of a numeric
quantity. It first gates on `is_numeric_quantity`: a non-numeric argument (a bare
symbol, `Positive[x]`) returns `NULL`, so the call is left unevaluated and the
symbolic expression flows on. For a numeric argument the shared helper
`numeric_real_sign` classifies reality and sign — exact `EXPR_INTEGER` /
`EXPR_BIGINT` / `Rational` / `EXPR_REAL` (and `EXPR_MPFR`) are read directly with
`expr_numeric_sign`, everything else is `numericalize`d at machine precision
(`numeric_machine_spec()`). A value that numericalizes to a genuinely non-real
complex number answers `False`; otherwise the verdict is `sign > 0`.

**Data structures.** Everything is `Expr`. A packed list or `NDArray` takes a
one-pass fast path: `ndint_sign_predicate(res, NDSP_POSITIVE)` (`src/ndinteger.c`)
reads the machine buffer and emits a `List` of `True`/`False` with no per-element
`Expr` and no evaluator round-trip. It declines the shapes it cannot cover —
complex buffers, rank > 1, and `Indeterminate` elements (unordered under IEEE) —
and `ndarray_delist_and_reeval` then threads the ordinary scalar path, which is
reached anyway by the `Listable` attribute for an unpacked list.

**Complexity / limits.** `O(1)` per scalar, `O(n)` element-wise over a buffer.
The verdict for an inexact or symbolic-constant argument rests on a
machine-precision numericalization, so a quantity indistinguishable from zero at
machine precision is classified by its machine sign. Attributes `Listable`,
`Protected`.
