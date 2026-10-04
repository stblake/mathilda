---
source: src/special_functions/sinhintegral.c
references:
  - "DLMF §6.2 — the hyperbolic sine integral Shi."
---
**Algorithm.** `builtin_sinhintegral` evaluates `Shi(z) = Int_0^z Sinh[t]/t dt`,
entire and odd, the imaginary-axis sibling of `Si` (`Shi(z) = -i Si(i z)`). Exact
special values: `0 -> 0`, `±Infinity`, `±I Pi/2` (directed imaginary infinity),
`ComplexInfinity`/`Indeterminate -> Indeterminate`; a negative-leading `Times`
folds by odd symmetry. Numeric (MPFR): the convergent Maclaurin series — whose
real-axis terms are all positive, so there is **no** catastrophic cancellation
and only a small fixed guard is needed (complex inputs keep the `~ |z|/ln2`
guard); and for large `|z|` the asymptotic `cosh(z) F(z) + sinh(z) G(z)` summed
to its smallest term, with the Stokes constant `i (Pi/2) sign(Im z)` restored
off the real axis. A machine-real result that overflows a `double`
(`Shi` grows like `e^|x|`) is emitted as a 53-bit MPFR real. The machine kernel
`sinhintegral_machine_complex` carries a `sf_series_usable` cancellation gate.

**Data structures.** `Expr`; `mpfr_t` (real) and `ncpx` (`mpfr_t` re/im,
complex); a `double complex` fallback for `USE_MPFR=0` builds. ND: unary kernel
`NDKU_SinhIntegral = { sinhintegral_machine_complex, ndk_SinhIntegral_r, ... }`
(the real kernel is `sf_machine_shi`), registered `REG_U`, so `packed_aware`.
Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Convergent regime needs no cancellation guard on the
real axis; off-axis and asymptotic regimes behave like `Si`.
`D[SinhIntegral[z], z] = Sinh[z]/z`. `Compile[]` lowers at both scalar and
rank-1 array shapes (`Compiled -> True`).
