---
source: src/numbertheory/primeomega.c
---
**Algorithm.** `builtin_primeomega` factors `|n|` and returns `Omega(n)`, the number of
prime factors counted with multiplicity — the sum of the exponents in the prime
factorisation (`primeomega_from_exps`). This is the quantity `LiouvilleLambda` forms before
taking `(-1)^Omega`. It takes one positional argument plus an optional `GaussianIntegers`
rule; under `GaussianIntegers -> True`, or for a non-real Gaussian `n`, the count runs over
the Gaussian prime factors. `Omega(1) = Omega(-1) = 0`.

**Data structures.** Factorisation via `df_factor_mpz` / `df_gaussian_prime_factor` into
`mpz_t` / `unsigned long` arrays; the exponent sum gives a machine `Integer`. `PrimeOmega`
has no dedicated ND/packed/`Compile` kernel; `Listable` threading is handled by the
evaluator.

**Complexity / limits.** Dominated by factoring `|n|` through `FactorInteger`;
`df_factor_mpz` verifies each base's primality with 40 Miller–Rabin rounds and declines —
leaving the call unevaluated — rather than trust a composite cofactor. `n == 0` is left
unevaluated, and a wrong argument count emits `PrimeOmega::argt`.
