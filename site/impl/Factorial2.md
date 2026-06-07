---
source: src/numbertheory.c
---
`builtin_factorial2` computes the double factorial `n!! = n(n-2)(n-4)…` by an explicit step-2 product loop. For small non-negative `EXPR_INTEGER` (n ≤ 30) it accumulates in an `int64_t`; for larger integers and `EXPR_BIGINT` it accumulates in a GMP `mpz_t` and returns an `EXPR_BIGINT`. Special cases: `(-1)!! = 0!! = 1`; negative integers return `ComplexInfinity` (poles of the analytic continuation). A `BigInt` argument too large for `mpz_fits_ulong_p` returns NULL (left symbolic) rather than attempting an unbounded loop. Non-integer arguments return NULL so `Factorial2[x]` stays symbolic. The function does not use the Gamma-based continuation for non-integers.
