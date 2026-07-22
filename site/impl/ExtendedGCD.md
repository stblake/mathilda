---
references:
  - "D. E. Knuth, *The Art of Computer Programming, Vol. 2: Seminumerical Algorithms*, 3rd ed. (Addison-Wesley, 1997), §4.5.2."
source: src/numbertheory.c
---
**Algorithm.** `builtin_extendedgcd` computes a multi-argument Bézout identity, returning `{g, {r1, ..., rN}}` with `g = GCD[n1, ...]` and `g = sum r_i n_i`. It folds GMP's `mpz_gcdext` pairwise: at each step `gcd(running_g, a_i) = s*running_g + t*a_i`, every previously accumulated cofactor is scaled by `s` and the new cofactor `t` is appended. `running_g` starts at `0`, so the first step yields `|a0|` with cofactor `sign(a0)`, and GMP normalises `g` non-negative throughout — matching the conventional sign. `ExtendedGCD[]` is `{0, {}}`.

**Data structures.** A heap array of `count` `mpz_t` cofactor accumulators plus scalar `mpz_t` registers for `gcdext`; results pass through `expr_bigint_normalize`. Integer-only — machine ints auto-promote via `expr_to_mpz`. Inexact (Real/MPFR) arguments emit `ExtendedGCD::exact`; exact-but-non-integer rationals emit `ExtendedGCD::egcd`; both return `NULL`.

**Complexity / limits.** Dominated by the `gcdext` chain, `O(N · M(d) log d)` for `d`-digit inputs (`M` = GMP's multiplication cost).
