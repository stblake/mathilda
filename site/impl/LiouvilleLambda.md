---
references:
  - "G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — the Liouville function lambda and the prime-factor count Omega (Chapter 17)."
source: src/numbertheory/liouvillelambda.c
---
**Algorithm.** `builtin_liouvillelambda` factors `|n|` and returns
`lambda(n) = (-1)^Omega(n)`, where `Omega(n)` is the sum of the prime-factor exponents
(`liouville_from_exps`). It takes one positional argument plus an optional `GaussianIntegers`
rule; under `GaussianIntegers -> True`, or for a non-real Gaussian `n`, the count runs over
the Gaussian prime factorisation. lambda is completely multiplicative, with
`lambda(-n) = lambda(n)` (sign ignored) and `lambda(1) = 1` (the empty product).

**Data structures.** The factorisation is a `mpz_t*` / `unsigned long*` pair from
`df_factor_mpz` (ordinary) or `df_gaussian_prime_factor` (Gaussian); only the exponent sum
is needed, and the result is a machine `Integer` of `+1` or `-1`. Unlike its sibling
`MoebiusMu`, `LiouvilleLambda` has no dedicated ND/packed/`Compile` kernel; `Listable`
threading is handled by the evaluator, so `LiouvilleLambda[Range[n]]` maps element-wise over
the boxed list.

**Complexity / limits.** Dominated by factoring `|n|` through `FactorInteger` (trial
division + Pollard rho + ECM). `df_factor_mpz` confirms each base prime with 40 Miller–Rabin
rounds and declines — leaving the call unevaluated — rather than return a confident wrong
value on an unfactored composite cofactor (the failure that once made these counting
functions wrong on an 82-digit semiprime). `n == 0` is left unevaluated, and a wrong
argument count emits `LiouvilleLambda::argt`.
