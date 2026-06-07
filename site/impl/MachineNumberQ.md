---
source: src/numeric.c
---
`builtin_machinenumberq` (`src/numeric.c`) returns `True` when the argument is a finite `EXPR_REAL` (via `is_machine_real_leaf`, which checks `EXPR_REAL` and `isfinite`), or a `Complex` whose real and imaginary parts are both finite machine reals; otherwise `False`. Exact integers/rationals and arbitrary-precision `EXPR_MPFR` values are not machine numbers.
