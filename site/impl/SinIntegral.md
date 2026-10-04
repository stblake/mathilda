---
source: src/special_functions/sinintegral.c
references:
  - "DLMF §6.2 — the sine integral Si."
  - "DLMF §6.6.5 — convergent Maclaurin series; §6.12.3 — asymptotic expansion."
---
**Algorithm.** `builtin_sinintegral` evaluates `Si(z) = Int_0^z Sin[t]/t dt`,
entire and odd. Exact special values: `0 -> 0`, `±Infinity -> ±Pi/2`,
`±I Infinity` (directed), `ComplexInfinity`/`Indeterminate -> Indeterminate`; a
negative-leading `Times` folds by odd symmetry. Numeric (MPFR): the convergent
Maclaurin series (DLMF 6.6.5) for small `|z|`, with `~ |z|/ln2` guard bits added
to absorb the cancellation exactly (partial sums reach `~ e^|z|`); the
asymptotic expansion (DLMF 6.12.3) `Si(z) = Pi/2 - cos(z) f(z) - sin(z) g(z)`
summed to its smallest term for large `|z|`. Complex arguments use the same two
regimes in the `ncpx` toolkit. The shared machine kernel
`sinintegral_machine_complex` carries a `sf_series_usable` cancellation gate
(declining to the MPFR implementation rather than returning a lossy answer).

**Data structures.** `Expr`; `mpfr_t` (real) and `ncpx` (`mpfr_t` re/im,
complex); a `double complex` fallback for `USE_MPFR=0` builds. ND: unary kernel
`NDKU_SinIntegral = { sinintegral_machine_complex, ndk_SinIntegral_r, ... }` (the
real kernel is `sf_machine_si`), registered `REG_U`, so `packed_aware`.
Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Convergent regime adds `~ |z|/ln2` guard bits;
asymptotic regime uses optimal truncation. Entire (defined for all `z`, real and
complex). `Compile[]` lowers at both scalar and rank-1 array shapes
(`Compiled -> True`).
