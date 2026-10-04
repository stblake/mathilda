---
references:
  - "T. M. Apostol, *Introduction to Analytic Number Theory*, Springer, 1976 — multiplicative functions and the divisor function sigma_k (Chapter 2)."
source: src/numbertheory/divisorsigma.c
---
**Algorithm.** `builtin_divisorsigma` takes `k` and `n` (plus an optional
`GaussianIntegers` rule), factors `|n|` via `df_factor_mpz` (which delegates to
`FactorInteger`), and builds the multiplicative formula
`sigma_k(n) = prod_i (p_i^((e_i+1) k) - 1) / (p_i^k - 1)` as an `Expr` tree
(`ds_build_factor`) that it then evaluates. One path serves integer, rational, radical and
fully symbolic `k` — `DivisorSigma[k, 12]` returns the product in closed form in `k`. The
special case `k == 0` degenerates to the divisor count `prod_i (e_i + 1)`
(`ds_divisor_count`). With `GaussianIntegers -> True`, or a non-real Gaussian `n`, it
factors into Gaussian primes (`df_gaussian_prime_factor`), normalises each prime to its
first-quadrant associate, and runs the same product.

**Data structures.** Primes and exponents come back as parallel `mpz_t*` / `unsigned long*`
arrays; prime atoms become `Expr` (`Integer`/`BigInt`, or `Complex[u, v]` for a Gaussian
prime), and the multiplicative factors are assembled as `Power`/`Plus`/`Times` trees and
reduced through `eval_and_free`. An exact-integer ND kernel exists: `ndk_DivisorSigma_ii`
in `src/ndinteger.c` (registered with `symtab_set_ndarray_binary_kernel`, so `DivisorSigma`
sits on both the `AWARE` and `INT64_OK` lists in `src/pack.c`) evaluates
`prod_i (1 + p^k + ... + p^(e_i k))` entirely in `int64` for a packed or visible `int64`
array, declining `n == 0`, `k < 0` (sigma_-1 is a `Rational`, which no int64 buffer holds)
and any overflow back to the GMP List path. `Compile[]` lowers `DivisorSigma[k, v]` at a
rank-1 integer-array shape (`Compiled -> True`) but not at a scalar shape — the
integer-only kernel has no scalar opcode.

**Complexity / limits.** Dominated by factoring `|n|` (trial division + Pollard rho + ECM
through `FactorInteger`); `df_factor_mpz` confirms every base prime with 40 Miller–Rabin
rounds and declines — leaving the call unevaluated — rather than trust a composite cofactor
on a hard semiprime. Sign of `n` is ignored, `n == 0` is left unevaluated, and a wrong
argument count emits `DivisorSigma::argrx`.
