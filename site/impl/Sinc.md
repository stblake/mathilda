---
source: src/special_functions/sinc.c
---
**Algorithm.** `builtin_sinc` evaluates the cardinal sine `Sinc[z] = Sin[z]/z`
with the removable singularity `Sinc[0] = 1`; it is entire and even. Each
argument kind takes the cheapest route: exact special values (`0 -> 1`,
`±Infinity -> 0`, `ComplexInfinity`/`Indeterminate -> Indeterminate`); a machine
real through libm (or `mpfr_sin(x)/x`); an arbitrary real through `mpfr_sin` at
the input precision; a complex argument through `sin(z)/z` in the shared `ncpx`
MPFR-complex toolkit; everything else stays symbolic.

**Data structures.** `Expr`; `mpfr_t` (real) and `ncpx` (`mpfr_t` re/im,
complex); a `double`/`double complex` fallback for `USE_MPFR=0` builds. ND:
real-only unary kernel `NDKU_Sinc = { NULL, ndk_Sinc_r, ... }` (the real kernel
is `sf_machine_sinc`), registered `REG_U`, so `packed_aware`. Attributes:
`Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** `O(1)` per element; the origin is handled exactly as
`1`. Symbolic for non-numeric arguments. `Compile[]` lowers at both scalar and
rank-1 array shapes (`Compiled -> True`).
