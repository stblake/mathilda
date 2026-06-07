---
source: src/precision.c
---
`builtin_set_precision` is a two-argument wrapper: it parses the precision argument into a `NumericSpec` via `parse_prec_arg` (accepting an integer/real digit count or `MachinePrecision`) and drives `numericalize(value, spec)` — the same engine `N` uses — to re-represent the value to the requested number of significant digits (an `EXPR_MPFR` at the corresponding bit width when MPFR is built, otherwise machine `double`). Returns `NULL` if the precision argument is not a valid spec.
