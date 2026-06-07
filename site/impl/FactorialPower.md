---
source: src/numbertheory.c
---
`builtin_factorialpower` computes the falling factorial `n^(k) = n(n-1)...(n-k+1)`. It requires `k` to be a concrete non-negative integer (`k == 0` gives `1`; negative or symbolic `k` returns `NULL`). For integer/bignum `n` it accumulates the exact product in GMP (`mpz_mul` of `n-i`), normalising via `expr_bigint_normalize`. For symbolic `n` with small `k` (`<= 32`) it expands the literal product `Times[n, n-1, ..., n-k+1]` through `eval_and_free` so `Expand`/`D` can act on it. The 1-argument step-form is not handled here.
