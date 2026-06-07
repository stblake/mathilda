---
source: src/complex.c
---
`builtin_sign` returns the sign (-1/0/1) of a real number — direct comparisons for `EXPR_INTEGER`/`EXPR_REAL` (and Rational by sign of numerator×denominator), `mpz_sgn` for BigInt, `mpfr_sgn` for MPFR. For a numeric `Complex[re, im]` with both parts numeric it returns the unit-modulus direction `z/Abs[z]` (short-circuiting `0+0I -> 0`); MPFR components take a fast path computing the direction directly via `mpfr_hypot` and division at the combined working precision rather than building the symbolic `z·Power[Abs[z], -1]` tree. Non-numeric arguments return `NULL` (unevaluated).
