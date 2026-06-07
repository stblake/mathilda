---
source: src/int.c
---
`builtin_integerexponent` returns the largest `k` with `b^k | n` (the base-`b` valuation, default base 10 = trailing-zero count). For base 2 it uses `mpz_scan1` (position of the lowest set bit = 2-adic valuation); otherwise GMP's `mpz_remove(q, |n|, base)`, which divides out `base` repeatedly and returns the count in one library call (`intexp_count`). `IntegerExponent[0, b]` is `Infinity` (every power divides 0). Validates arity (`::argt`), numeric non-integer `n` (`::int`), and base `>= 2` (`::ibase`); symbolic `n` returns `NULL`.
