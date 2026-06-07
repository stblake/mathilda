---
source: src/core.c
---
`builtin_oddq` (`src/core.c`) returns `True` for an `EXPR_INTEGER` with `n % 2 != 0`, uses `mpz_odd_p` for an `EXPR_BIGINT`, and returns `False` for everything else.
