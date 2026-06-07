---
source: src/numbertheory.c
---
**Algorithm.** `builtin_gcd` folds the arguments pairwise. It classifies them in one pass and chooses a path: an `int64` fast path using the binary/Euclidean `gcd`/`lcm` helpers; a GMP path (`mpz_gcd`) when any argument is a `EXPR_BIGINT`; and a rational fold for rational-like inputs using the identity `gcd(a/b, c/d) = gcd(a,c)/lcm(b,d)`, accumulating numerator with `mpz_gcd` and denominator with `mpz_lcm`. `GCD[]` is `0`, `GCD[x]` is `|x|`. All numerators/denominators are taken in absolute value before folding; any non-rational argument makes the call return `NULL` (left symbolic).

**Data structures.** Pure GMP `mpz_t` running accumulators; results pass through `expr_bigint_normalize` to demote back to `EXPR_INTEGER` when they fit, and `mpz_pair_to_rational_expr` reduces a num/den pair (dividing by their `mpz_gcd`) into an `Integer` or canonical `Rational`. GMP's `mpz_gcd` uses a subquadratic (HGCD) algorithm.
