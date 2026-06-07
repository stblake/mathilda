---
source: src/core.c
---
`builtin_evenq` (`src/core.c`) returns `True` for an `EXPR_INTEGER` with `n % 2 == 0`, uses `mpz_even_p` for an `EXPR_BIGINT`, and returns `False` for everything else.
