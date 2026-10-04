---
references:
  - "G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — the Moebius function and Moebius inversion (Chapter 16)."
source: src/numbertheory/moebiusmu.c
---
**Algorithm.** `builtin_moebiusmu` takes exactly one argument, factors `|n|`, and returns
`mu(n)`: `0` if any exponent is `>= 2` (a squared prime factor), otherwise `(-1)^m` for `m`
distinct primes, with `mu(1) = 1` (`moebiusmu_from_exps`). A non-real Gaussian-integer
argument is auto-detected and factored over `Z[i]` (the unit factor does not count). The
sign of `n` is ignored, matching `mu(-n) = mu(n)`.

**Data structures.** Factorisation via `df_factor_mpz` / `df_gaussian_prime_factor` into
`mpz_t` / `unsigned long` arrays; the answer is a machine `Integer` in `{-1, 0, 1}`. An
exact-integer ND kernel exists: `ndk_MoebiusMu_ii` in `src/ndinteger.c` (registered with
`symtab_set_ndarray_unary_kernel`, placing `MoebiusMu` on both the `AWARE` and `INT64_OK`
lists in `src/pack.c`) factors each `int64` element in place and declines `n == 0` or a
factor past the trial-division ceiling back to the GMP List path. `Compile[]` lowers
`MoebiusMu[v]` at a rank-1 integer-array shape (`Compiled -> True`) but not at a scalar
shape — the integer-only kernel is a real path over an array and no path at all over a
scalar (the `narrowing_only` branch in `src/compile/compile.c`).

**Complexity / limits.** Dominated by factoring `|n|`. As with its siblings, `df_factor_mpz`
verifies each base's primality with 40 Miller–Rabin rounds and declines rather than trust a
composite cofactor — which is what stopped `MoebiusMu` of an 82-digit semiprime returning a
confident wrong `-1`. `MoebiusMu[0]` is left unevaluated, and a wrong argument count emits
`MoebiusMu::argx`.
