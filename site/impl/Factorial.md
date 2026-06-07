---
source: src/numbertheory.c
---
**Algorithm.** `builtin_factorial` reduces only concrete numeric arguments. A non-negative integer `n` uses an `int64` loop for `n <= 20` and GMP's `mpz_fac_ui` beyond that; a negative integer gives `ComplexInfinity` (pole of `Gamma`). A machine `Real` evaluates `tgamma(x+1)`; an MPFR real evaluates `mpfr_gamma` of `x+1` at the input precision. Half-integer arguments (`d == ±2`) are folded to the closed form `coeff * Sqrt[Pi]` via the double-factorial relation, building `coeff` as an exact `Rational`. Other rationals and symbolic inputs return `NULL`. A `EXPR_BIGINT` argument is deliberately left symbolic — its factorial is astronomically large.

**Data structures.** GMP `mpz_t` for the bignum branch; the half-integer branch assembles `Times[Rational[...], Power[Pi, 1/2]]` and reduces it through `eval_and_free`.
