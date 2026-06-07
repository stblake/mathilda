---
source: src/core.c
---
`builtin_atomq` (`src/core.c`) returns `True` for any non-`EXPR_FUNCTION` node (integers, reals, bigints, symbols, strings) and for the two function heads Mathilda treats as atomic, `Rational` and `Complex`; every other `EXPR_FUNCTION` yields `False`.
