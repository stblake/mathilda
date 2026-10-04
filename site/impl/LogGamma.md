---
source: src/special_functions/loggamma.c
references:
  - "DLMF §5.11 — Stirling's series for ln Gamma."
  - "C. Lanczos, A precision approximation of the gamma function, SIAM J. Numer. Anal. 1 (1964) 86-96."
---
**Algorithm.** `builtin_loggamma` evaluates the analytic continuation of
`log Gamma(z)` (the continuous branch, **not** `Log[Gamma[z]]`). Exact: a
positive integer `n` gives `Log[(n-1)!]`; a half-integer gives the log of the
exact `Sqrt[Pi]` form plus the branch term `-Ceiling[-z] Pi I` for `z < 0`; a
non-positive integer is the pole `Infinity`; the symbolic infinities map to
`Infinity`/`ComplexInfinity`/`Indeterminate`. Numeric: machine real `z > 0` uses
libm `lgamma`; `z < 0` non-integer uses `lgamma` plus the imaginary branch term
`Im = -Pi Ceiling[-z]` (a complex result); arbitrary real uses `mpfr_lgamma`
plus the branch term; machine complex uses a `(g = 7, n = 9)` Lanczos log-gamma
(with a winding-correct reflection for `Re z < 1/2`); arbitrary complex uses the
Stirling asymptotic series with integer argument reduction and exact Bernoulli
numbers. Intervals route through `interval_apply_function`.

**Data structures.** `Expr`; `double complex` Lanczos; a local `lcx` (`mpfr_t`
re/im pair) toolkit for the Stirling path with an exact `mpq` Bernoulli
generator; GMP. ND: unary kernel `NDKU_LogGamma = { loggamma_machine_complex,
ndk_LogGamma_r, ... }` — the real part via `lgamma` for `x > 0`, the complex part
the *continued* log-gamma — registered `REG_U`, so `packed_aware`. Attributes:
`Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Lanczos is `O(1)`; the Stirling path uses `K ~ wp/2`
Bernoulli terms after shifting `Re(w)` up to `~ K/3`, so cost grows with the
requested precision. The imaginary part grows without bound with `Im z` (the
whole reason this is not `Log[Gamma[z]]`). `Compile[]` lowers at both scalar and
rank-1 array shapes (`Compiled -> True`).
