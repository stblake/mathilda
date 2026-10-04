---
source: src/special_functions/polygamma.c
references:
  - "DLMF §5.15 — the polygamma functions."
  - "M. Abramowitz & I. A. Stegun, Handbook of Mathematical Functions, §6.4.11 — the asymptotic expansion."
---
**Algorithm.** `builtin_polygamma` evaluates `PolyGamma[n, z] = psi^(n)(z)`;
`PolyGamma[z]` rewrites to `PolyGamma[0, z]`. Order `-1` emits the inert
`LogGamma[z]`; `z` a non-positive integer is the pole `ComplexInfinity`. At a
positive integer `z`: `n = 0` gives `H_{z-1} - EulerGamma`; odd `n >= 1` gives a
rational plus rational·`Pi^(n+1)` (via exact Bernoulli / zeta); even `n >= 2`
stays symbolic (zeta of an odd argument). At a positive rational the digamma
(`n = 0`) uses the Gauss digamma theorem — but only half-integers auto-evaluate,
general rationals being routed to the numeric path to avoid catastrophic
cancellation. Numeric: an inexact real uses `mpfr_digamma` for `n = 0`, else the
recurrence-shift + Bernoulli asymptotic series (A&S 6.4.11) in `mpfr`; an inexact
complex argument runs the same asymptotic in complex arithmetic. Intervals route
through `interval_polygamma`.

**Data structures.** `Expr`; an exact `mpq` Bernoulli cache (process-lifetime);
a local `pcx` (`mpfr_t` re/im) toolkit for the asymptotic; GMP. ND: unary kernel
`NDKU_PolyGamma` (real digamma via `sf_machine_digamma`) **and** binary kernel
`NDKB_PolyGamma = { ndk_PolyGamma_c, ... }` (order and argument), registered
`REG_U` and `REG_B`, so `packed_aware` — the evaluator canonicalises
`PolyGamma[x]` to `PolyGamma[0, x]`, so the binary kernel carries the real work.
Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Exact-argument cap `100000`, exact-order cap `256`,
numeric-order cap `1024`. The asymptotic series needs the argument shifted up to
`Re ~ 0.13 wp + n`; cost grows with the requested precision. `Compile[]` lowers
both `PolyGamma[x]` and `PolyGamma[n, x]` at scalar and rank-1 array shapes
(`Compiled -> True`).
