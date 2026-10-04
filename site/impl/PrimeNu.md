---
source: src/numbertheory/primenu.c
---
**Algorithm.** `builtin_primenu` factors `|n|` and returns `nu(n)`, the number of distinct
prime factors — the count of factors, independent of their exponents (`primenu_from_count`).
It is the additive companion to `PrimeOmega` (which sums the exponents); `nu` and `Omega`
coincide exactly when `n` is square-free. It takes one positional argument plus an optional
`GaussianIntegers` rule; under `GaussianIntegers -> True`, or for a non-real Gaussian `n`,
the count runs over the distinct Gaussian prime factors. `nu(1) = nu(-1) = 0`.

**Data structures.** Factorisation via `df_factor_mpz` / `df_gaussian_prime_factor` into
`mpz_t` / `unsigned long` arrays; only the factor count is used, and the result is a machine
`Integer`. Unlike `MoebiusMu`, `PrimeNu` has no dedicated ND/packed/`Compile` kernel;
`Listable` threading is handled by the evaluator.

**Complexity / limits.** Dominated by factoring `|n|` through `FactorInteger`;
`df_factor_mpz` confirms each base prime with 40 Miller–Rabin rounds and declines — leaving
the call unevaluated — on an unfactored composite cofactor rather than returning a wrong
count. `n == 0` is left unevaluated, and a wrong argument count emits `PrimeNu::argt`.
