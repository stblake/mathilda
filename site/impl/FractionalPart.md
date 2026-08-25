---
source: src/piecewise.c
---
`builtin_fractionalpart` computes `x - IntegerPart[x]` through the shared `do_piecewise(res, OP_FRACPART, ...)` kernel, keeping the sign of `x` and the precision of the input: `EXPR_REAL` returns `v - trunc(v)`, `EXPR_MPFR` subtracts `mpfr_trunc` at full precision, and exact rationals return an exact `Rational`. Registered `PROTECTED | NUMERICFUNCTION | LISTABLE`; a quantity with no monotone reduction (e.g. `FractionalPart[10^7 3^(2/3)]`) is left symbolic.
