# Task: Eliminate inverse-function substitution pass

## Goal
`Eliminate[{Dt[y]==(x ArcSin[x])/Sqrt[1-x^2] Dt[x], u==ArcSin[x], Dt[u]==Dt[x]/Sqrt[1-x^2]}, {x, Dt[x]}]`
must return `u Dt[u] Sin[u] == Dt[y]` (FullForm `Equal[Times[u, Dt[u], Sin[u]], Dt[y]]`).
Fully general across all inverse transcendental functions.

## Root cause (verified)
`ArcSin[x]` buried inside eq1 is an opaque function-atom containing elim var `x`
=> `nlin` bail. `try_inverse_rewrite` only inverts a top-level `f[u]==v`; it never
propagates `ArcSin[x] -> u` into other equations.

## Design (validated empirically)
Feeding the hand-transformed system
`{Dt[y]==x u Dt[x]/Cos[u], Dt[u]==Dt[x]/Cos[u], x==Sin[u]}` into current Eliminate
already yields the exact target. Add a pre-pass that performs that transform.

For a defining eq `M == ArcF[x]` (M elim-free, x single elim symbol, ArcF invertible):
1. globally replace `ArcF[x] -> M`;
2. globally replace companion radical `(comp_base)^(p/2) -> coFn[M]^(co_sign*p)`
   (keeps denominators as main-var monomials -> gb_poly_strip_monomial divides out);
3. replace defining-eq slot with `x == F[M]`;
4. emit Eliminate::ifun.

### Identity table
ArcSin  Sin   1-x^2   Cos   +1
ArcCos  Cos   1-x^2   Sin   +1
ArcTan  Tan   1+x^2   Cos   -1  (Sec)
ArcSinh Sinh  1+x^2   Cosh  +1
ArcCosh Cosh  x^2-1   Sinh  +1
ArcTanh Tanh  1-x^2   Cosh  -1  (Sech)
Log     Exp   none    none      (main-var E^M -> nlin downstream, known limit)

## Steps
- [ ] Implement pre-pass in src/poly/eliminate.c
- [ ] Wire in before forward transcendental pass
- [ ] Verify headline example exact target
- [ ] No regression in eliminate_tests
- [ ] Add extensive unit tests
- [ ] Update docs/spec + changelog
- [ ] valgrind smoke

## Review
(to fill in)

## Review (DONE 2026-07-11)
Implemented `inv_substitute_all` + helpers in src/poly/eliminate.c, wired in
before the forward transcendental pass. All 6 inverse families + both Dt forms
+ shorthand verified. Exact target `Equal[Times[u, Dt[u], Sin[u]], Dt[y]]`.
- eliminate_tests: PASS (13 new tests, 0 regressions)
- integrate_derivdivides_tests: PASS (new test_inverse_trig_substitution)
- integrate_dispatch_tests: PASS
- Unlocked Integrate[x ArcSin[x]/Sqrt[1-x^2], x] = x - Sqrt[1-x^2] ArcSin[x]
  (+ ArcCos/ArcSinh/ArcTan kin), D[]-verify to 0.
- leaks-clean (0 bytes / 30 iters); strict C99 -Wall -Wextra clean.
- Docs: structural-manipulation.md bullet + changelog 2026-07-06.md.

### Correction handled mid-task
Initial table included Log -> regressed test_log_power (`u==Log[x]` was
intercepted and substituted x->E^M, an un-atomizable main-var exponential).
Fix: exclude Log; the forward exp/log pass already handles it.
