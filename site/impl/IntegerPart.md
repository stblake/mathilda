---
source: src/piecewise.c
---
`builtin_integerpart` dispatches through the shared `do_piecewise(res, OP_INTPART, ...)` kernel, truncating each numeric type *toward zero*: Integer/BigInt pass through; rationals use `mpz_tdiv_q` on numerator and denominator; `EXPR_REAL` uses C `trunc()`; `EXPR_MPFR` uses `mpfr_trunc` then `mpfr_get_z` into an `mpz_t`, so arbitrarily large values never silently overflow `int64_t`. Registered `PROTECTED | NUMERICFUNCTION | LISTABLE`. It satisfies `IntegerPart[x] + FractionalPart[x] == x`; non-numeric arguments stay symbolic.
