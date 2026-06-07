---
source: src/core.c
---
`builtin_numberq` (`src/core.c`) returns `True` for an explicit number — `EXPR_INTEGER`, `EXPR_REAL`, `EXPR_BIGINT`, `EXPR_MPFR` (under `USE_MPFR`), or a `Rational`/`Complex` head — and `False` otherwise. (Contrast `NumericQ`, whose `is_numeric_quantity` helper also accepts symbolic constants like `Pi` and numeric-function calls.)
