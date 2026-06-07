---
source: src/numbertheory.c
---
**Algorithm.** `builtin_lcm` mirrors `builtin_gcd`. It folds pairwise with `lcm(a,b)=ab/gcd(a,b)`: an `int64` fast path, a GMP path (`mpz_lcm`) when any argument is a `EXPR_BIGINT`, and a rational fold using `lcm(a/b, c/d) = lcm(a,c)/gcd(b,d)` (numerator accumulated with `mpz_lcm`, denominator with `mpz_gcd`). A zero argument zeroes the running LCM (and short-circuits). `LCM[]` is `1`, `LCM[x]` is `|x|`; non-rational arguments return `NULL`.

**Data structures.** GMP `mpz_t` accumulators; `expr_bigint_normalize` demotes results that fit in `int64`, and `mpz_pair_to_rational_expr` reduces the rational result. Shares the rational num/den coercion helpers with GCD.
