---
source: src/facint.c
---
`builtin_nextprime` computes `NextPrime[x]` / `NextPrime[x, k]`. The start point is taken as `⌊x⌋` for Real/Rational and exactly for Integer/BigInt, into a GMP `mpz_t`. `k = 0` returns x unchanged. For `k > 0` it iterates GMP's `mpz_nextprime` k times (the next probable prime strictly greater than the current value). For `k < 0` it iterates `mpz_prevprime` |k| times, returning unevaluated (NULL) if it would step at or below 2 or no previous prime exists. The result is normalised via `expr_bigint_normalize`.
