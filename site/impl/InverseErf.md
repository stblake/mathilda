---
source: src/special_functions/inverf.c
references:
  - "DLMF §7.17 — inverse error functions."
  - "S. Winitzki, A handy approximation for the error function and its inverse (2008) — the closed-form seed."
---
**Algorithm.** `builtin_inverf` solves `s = erf(z)`. Per Mathematica, explicit
values are produced only for *real* `s` in `[-1, 1]`; complex and `|s| > 1`
inputs stay symbolic. Exact: `0 -> 0`, `1 -> Infinity`, `-1 -> -Infinity`; a
symbolic negative-leading `Times` folds by odd symmetry `InverseErf[-x] ->
-InverseErf[x]`. Machine real `|s| < 1` uses `inverf_double`: a Winitzki
closed-form seed polished by four Newton steps on the libm `erf`,
`z <- z - (erf(z) - s)(sqrt(pi)/2) e^{z^2}`. Arbitrary-precision real uses Newton
with precision doubling on `mpfr_erf`, the final full-precision `erf` dominating
cost. The two-argument `InverseErf[z0, s]` rewrites to
`InverseErf[s + Erf[z0]]`, staying symbolic if `Erf[z0]` does not reduce. The
derivative `(sqrt(pi)/2) e^{InverseErf[z]^2}` lives in `calculus/deriv.c`; Series
follows via Taylor-via-`D`.

**Data structures.** `Expr`; `double` (`inverf_double`, exposed for ND and the
`Compile[]` VM); `mpfr_t` for the arbitrary-precision Newton loop. ND: real-only
unary kernel `NDKU_InverseErf = { NULL, ndk_InverseErf_r, ... }` (no complex
kernel), registered `REG_U`, so `packed_aware`. Attributes: `Listable`,
`NumericFunction`, `Protected`.

**Complexity / limits.** Newton converges quadratically; the MPFR cost tracks
the final full-precision `erf`. Domain is the real interval `[-1, 1]`; outside
it (and for complex or exact-rational `s`) the call stays symbolic.
`Compile[]` lowers at both scalar and rank-1 array shapes (`Compiled -> True`).
