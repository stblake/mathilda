---
source: src/core.c
---
`builtin_length` (in `src/core.c`) returns an `EXPR_INTEGER` equal to the argument's `arg_count` when it is an `EXPR_FUNCTION`, and `0` for atoms (which have no parts).
