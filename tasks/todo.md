# TODO: DSolve M58 — higher-order autonomous reduction (missing-x, order n ≥ 3)

Plan: `/Users/user/.claude/plans/let-s-continue-our-implementation-binary-hartmanis.md`
Lift `AutonomousReduction` (`dsolve_autonomous.c`) from order-2-only to any order n≥2.
Missing-y at order 3+ already done (`lower_reduce`); this is the missing-x gap.
Predicted: all 7 Group-A reduced ODEs close (263/264/267/268/1143/1167/1168).

## Build
- [ ] Generalize `dsolve_autonomous_try`: gate `>=2`; derivative chain D[1..n]
      (D[k+1]=p·d/dy(D[k])); eq1 = D[n]==F[y^(k)→D_k]; `dsolve_renumber_constants`
      freeze C[1..n-1]→C[2..n]; stage-2 separable + degenerate rejection unchanged.
- [ ] Zero cascade/wiring changes (method already registered). Update docstring.
- [ ] `make -j` + `make check-c99` green.

## Verify
- [ ] ORDER-2 REGRESSION FIRST: `y y''==(y')^2` → `C[2]E^(C[1]x)`; existing autonomous tests.
- [ ] Order-3 spot-check: 1143/1167/1168/263/264/267/268 solve + back-substitute ~0.
- [ ] Units `t_m58_*` (test_dsolve.c, isolation). Stress `test_dsolve_m58_stress.c` + CMake.
- [ ] All DSolve stress suites (m5/m12/m14/m18/m55/m56/m57/m58) green.
- [ ] CLEAN corpus run (no concurrent Mathilda): §2.1.2 delta, 0 FAIL, lower baseline 629.
- [ ] valgrind spot-check one order-3 solve (ownership clean in the file).

## Docs + release
- [ ] calculus.md (AutonomousReduction → order n≥2), changelog 2026-09-21.md (M58),
      DSOLVE_PLAN.md (M58 block + §1d), STATUS.md + regenerate reports/2.1.2.{tsv,md}.
- [ ] version.h 0.173→0.174; commit `; v0.174`; tag v0.174.

## Review
M58 complete (v0.174). `DSolve`AutonomousReduction` lifted from order-2-only to any order n≥2: the
derivative chain `D_{k+1}=p·d/dy(D_k)` reduces `y⁽ⁿ⁾==f(y,…,y⁽ⁿ⁻¹⁾)` (missing x) to an order-(n−1)
ODE in p(y), then `y'==p` separable; stage-1 constants frozen C[1..n−1]→C[2..n] via
`dsolve_renumber_constants`. Zero cascade/wiring changes — confined to `dsolve_autonomous.c`. For n=2
the chain reproduces the classical `p·p_y==f` exactly (order-2 exp/Tan-Tanh solve, elliptic declines
— verified).

**The snag (surfaced to user, who chose "ship the safe increment"):** stage-1 (reduction) always
closes, but stage-2's separable quadrature `∫dy/p` is non-elementary for most order-3 cases, and
Integrate SPINS uninterruptibly there (TimeConstrained can't preempt — 1167/264 hit 45 s). Added a
stage-2 spin-guard (decline a Log/y-denominator radicand before the spin: 45 s → 0.1 s) + decline
memo + 5 s deadline. So §2.1.2 yield is **+1** (1143, `p=Sqrt[C₁y²+C₂]` elementary), not the plan's
+7 — stage-1 closure ≠ stage-2 closure. The full order-3 set needs an implicit-first-integral output
(inert quadrature) — deferred as the "big lever" future item.

§2.1.2 clean re-run **590 → 591, 0 FAIL** (deterministic +1 = 1143; net P↔U is the flaky
`_with_linear_symmetries` cluster). Gate baseline 629 → 628. Stress m5/m12/m56/m57 + new m58 green;
autonomous cases in dsolve_stress pass (undetcoeff abort is pre-existing); check-c99 green; 0 crashes
across 1204 corpus cases (memory clean). Order-2 constant labels are cosmetically swapped
(C[1]E^(C[2]x)) — tests use back-substitution + C[2]-presence, so unaffected.

Lessons: (1) for a reduction method, stage-1 solvability does NOT imply the follow-on quadrature
closes — must check the WHOLE pipeline before estimating yield. (2) TimeConstrained cannot preempt an
uninterruptible Integrate spin — gate BEFORE the spin, don't time-bound it. (Recorded in memory.)

Pre-existing (flagged, not M58): `dsolve_tests` t_m19 in-suite abort; `dsolve_stress` t_stress_undetcoeff
abort (both on clean HEAD). The M56 first-integral method now solves `y''==2y³` (implicit first integral),
so the old `t_auto_declines_elliptic` was stale — updated to assert the PINNED autonomous method declines.
