---
source: src/special_functions/gamma.c
references:
  - "DLMF §5.2, §5.5 — the Euler gamma function and the reflection formula."
  - "DLMF §8.2, §8.11 — the (upper) incomplete gamma function Gamma(a,z)."
  - "J. L. Spouge, *Computation of the gamma, digamma, and trigamma functions*, SIAM J. Numer. Anal. 31 (1994) 931–944."
  - "C. Lanczos, *A precision approximation of the gamma function*, SIAM J. Numer. Anal. Ser. B 1 (1964) 86–96."
---
**Algorithm.** `builtin_gamma` serves `Gamma[z]`, `Gamma[a, z]` (upper incomplete), and `Gamma[a, z0, z1] = Gamma[a,z0] - Gamma[a,z1]`. For `Gamma[z]` the layering is: an **exact integer / half-integer** reduces to `(z-1)!` through the `Factorial` machinery (exact `BigInt`, or the `Sqrt[Pi]` half-integer rationals); symbolic infinities map to `Infinity`/`ComplexInfinity`/`Indeterminate`; a **machine real** uses libm `tgamma` (overflow promotes to a 53-bit MPFR real via `mpfr_gamma`, poles give `ComplexInfinity`); an **arbitrary real** uses `mpfr_gamma`; a **machine complex** uses a fixed-coefficient **Lanczos** approximation (`g = 7, n = 9`, ~15 digits, reflection for `Re z < 1/2`); and an **arbitrary-precision complex** uses **Spouge's** approximation (runtime-computable coefficients, so it honours the requested precision) in the file-local complex-MPFR toolkit `gcx`. The incomplete form: `Gamma[a, 0] = Gamma[a]`, `Gamma[a, Infinity] = 0`; a **positive integer `a`** with symbolic/exact `z` expands to the finite closed form `(n-1)! e^{-z} Sum_{k<n} z^k/k!`; numeric real routes to `mpfr_gamma_inc`; numeric complex uses a lower-series / Lentz-continued-fraction split in `gcx` with the `e^{-z} z^a` prefactor.

**Data structures.** `Expr`; file-local `gcx` (pairs of `mpfr_t`, alias-safe, explicit precision) for complex/incomplete MPFR work; GMP for the exact factorial path. The ND kernel is a real `REG_U` registration (`NDKU_Gamma`): real buffers via `ndk_Gamma_r` (libm `tgamma`), complex buffers via `gamma_machine_complex` (the machine Lanczos kernel). `Compile[]` lowers `Gamma` at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes.

**Complexity / limits.** `O(1)` per element on the real/machine-complex paths; the MPFR paths scale their term count with precision. Arbitrary-precision complex uses Spouge (precision-honest); machine complex uses Lanczos (~15 digits). Simple poles at the non-positive integers give `ComplexInfinity`. The incomplete integer expansion is capped at `a <= 1000`. Attributes: `Listable`, `NumericFunction`, `Protected`.
