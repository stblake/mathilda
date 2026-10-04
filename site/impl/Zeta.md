---
source: src/special_functions/zeta.c
references:
  - "DLMF §25.2 — the Riemann zeta function; §25.11 — the Hurwitz zeta function."
  - "DLMF §25.11.5 — the Euler-Maclaurin summation formula used for the Hurwitz / complex kernel."
---
**Algorithm.** `builtin_zeta` evaluates `Zeta[s]` (Riemann) and `Zeta[s, a]`
(Hurwitz). Exact integer `s`: `1 -> ComplexInfinity` (pole), `0 -> -1/2`, even
`2n > 0 -> rational · Pi^(2n)` (via exact Bernoulli), negative `-m -> rational`
(Bernoulli; `0` at even `m`), odd `2n+1 > 0` left symbolic; `s = Infinity -> 1`.
Hurwitz exact: `a = 1 -> Zeta[s]`; `a = 1/2 -> (2^s - 1) Zeta[s]`; a positive
integer `a -> Zeta[s] - Sum_{k=1}^{a-1} k^-s`. Numeric: real Riemann zeta via
`mpfr_zeta`; complex `s`, or any `a != 1`, via a Euler-Maclaurin complex-MPFR
kernel (DLMF 25.11.5) — head sum plus the `(a+N)^(1-s)/(s-1)`, `(1/2)(a+N)^-s`
and Bernoulli correction terms, with `N` chosen from the working precision and
`|s|` and the correction series truncated at its optimal (smallest) term; the
two-argument form uses the *symmetric* power `((a+k)^2)^(-s/2)` for `Re a < 0`.
Intervals route through `interval_apply_function`.

**Data structures.** `Expr`; an exact `mpq` Bernoulli cache; a local `zcx`
(`mpfr_t` re/im) toolkit for the Euler-Maclaurin kernel; GMP; `mpfr_zeta` for
real Riemann. ND: unary kernel `NDKU_Zeta = { NULL, ndk_Zeta_r, ... }` (real
Riemann, via `sf_machine_zeta`), registered `REG_U`, so `packed_aware`. The
two-argument Hurwitz form has **no** ND buffer kernel on `Zeta` — the separate
`HurwitzZeta` head carries the binary kernel `NDKB_HurwitzZeta`. Attributes:
`Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Exact-integer cap `ZETA_EXACT_INT_CAP = 10000`;
Hurwitz-integer-`a` cap `100000`. The Euler-Maclaurin head term count is
`N ~ digits + |s|`, the correction series optimally truncated. Real zeta goes
only through `mpfr_zeta`. `Compile[]` lowers `Zeta[s]` at both scalar and rank-1
array shapes (`Compiled -> True`).
