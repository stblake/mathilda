---
source: src/int.c
---
`builtin_integerlength` returns the number of base-`b` digits of `|n|` (default base 10). For bases `<= 62` it uses GMP's `mpz_sizeinbase`, which is exact for power-of-two bases and at most one too large otherwise — corrected by comparing `|n|` against `base^(s-1)` (`intlen_count_digits`). For arbitrary-precision bases it counts via repeated `mpz_tdiv_q`. `IntegerLength[0]` is `0`. Validates arity (`IntegerLength::argt`), numeric non-integer `n` (`::int`), and base `>= 2` (`::ibase`); symbolic `n` returns `NULL`.
