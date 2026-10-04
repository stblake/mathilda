---
source: src/special_functions/inverfc.c
references:
  - "DLMF §7.17 — inverse error functions."
  - "S. Winitzki, A handy approximation for the error function and its inverse (2008) — the closed-form seed."
---
**Algorithm.** `builtin_inverfc` solves `s = erfc(z)`. Since `erfc` maps the real
line onto `(0, 2)`, explicit values are produced only for *real* `s` in `[0, 2]`;
out-of-domain and complex inputs stay symbolic. Exact: `0 -> Infinity`,
`1 -> 0`, `2 -> -Infinity`. Machine real `0 < s < 2` uses `inverfc_double`: a
Winitzki seed on `w = 1 - s` polished by four Newton steps **directly on the
libm `erfc`**, `z <- z + (erfc(z) - s)(sqrt(pi)/2) e^{z^2}` — it deliberately
does **not** route through `InverseErf[1 - s]`, which would lose all significance
to cancellation for small `s`. Arbitrary precision uses Newton with precision
doubling on `mpfr_erfc`. The derivative `-(sqrt(pi)/2) e^{InverseErfc[z]^2}`
lives in `calculus/deriv.c`; Series follows via Taylor-via-`D`. (`erfc` is not
odd, so there is no auto-applied symmetry rewrite.)

**Data structures.** `Expr`; `double` (`inverfc_double`, exposed for ND and the
`Compile[]` VM); `mpfr_t` for the arbitrary-precision Newton loop. ND: real-only
unary kernel `NDKU_InverseErfc = { NULL, ndk_InverseErfc_r, ... }`, registered
`REG_U`, so `packed_aware`. Attributes: `Listable`, `NumericFunction`,
`Protected`.

**Complexity / limits.** Newton converges quadratically; the MPFR cost tracks
the final full-precision `erfc`. Domain is the real interval `[0, 2]`; outside
it (and for complex or exact-rational `s`) the call stays symbolic.
`Compile[]` lowers at both scalar and rank-1 array shapes (`Compiled -> True`).
