# DSolve M20 — first-order symmetry gap, Stage 1: PolynomialShiftSubstitution

Plan: `/Users/user/.claude/plans/twinkly-chasing-hammock.md`
Method: `y' = R(x) + g(x)(φ(x)+c·y)^p` → `u=φ+c·y` → separable → implicit first integral.
Generalizes FirstOrderSubstitution (linear arg ax+by+c) to x-dependent polynomial shift.
Targets radical `[F(x),G(x)]` cases that abaco2_similar ABORTs: 365/371/372/376/378/402/403/424.

## Step 0 — study templates
- [ ] dsolve_fos.c (FirstOrderSubstitution: three-fn contract, substitution + recurse)
- [ ] dsolve_run_implicit substrate + dsolve_method_builtin registration
- [ ] cascade wiring in dsolve.c (enum, ds_method_from_string, extern, dispatch slot)

## Step 1 — dsolve_polyshift.c
- [ ] detect radical/power atom base^p (p non-integer), base linear in y → φ(x), c
- [ ] u=φ+c·y; G(x,u)=φ'+c·f(x,(u-φ)/c); require free of y (general substitution test)
- [ ] recurse dsolve_run on u'==G; back-substitute u=φ+c·y; emit implicit first integral
- [ ] node budget + transcendental guard → fast decline on non-match
- [ ] register DSolve`PolynomialShiftSubstitution + ATTR_PROTECTED + docstring; init hook

## Step 2 — cascade
- [ ] slot after FirstOrderSubstitution (fos), before LieSymmetry

## Step 3 — tests + gates + docs
- [ ] test_dsolve_m20_stress.c (φ,c,g,p forward grid) + CMake
- [ ] t_m20_* units in test_dsolve.c (402/371 solve + non-linear-base decline)
- [ ] all DSolve ctest suites + check-c99 green; valgrind decline-path clean
- [ ] re-measure §2.1.2: positive delta, 0 FAIL, 0 P→U regressions; lower gate baseline
- [ ] DSOLVE_PLAN.md M20 + §1a; calculus.md; changelog; STATUS.md

## Review

Landed **M20 Stage 1**: `DSolve\`PolynomialShiftSubstitution` (`dsolve_polyshift.c`) — the
x-dependent-shift generalisation of `FirstOrderSubstitution`. `y' = R(x) + g(x)(φ(x)+c y)^p`
(R=−φ'/c) → `u=φ+c y` → separable `u'=c g u^p` → implicit first integral (branch-safe).
Deterministic replacement for the radical `[F(x),G(x)]`-symmetry cases `abaco2_similar`
aborts on. Cascade slot: after Separable, before Linearizable/Exact/Lagrange.

**Solves 2.1.2-402/-371/-372/-376/-403/-424** (6 of the 8 radical sub-cluster);
365 (general separable) and 378 (x-dependent c) correctly decline.

Two robustness lessons (both caused hangs, both fixed):
1. `Simplify[(x+Sqrt[u])u^(-1/2)]` LOOPS → build the reduced form by plain evaluation, not
   `ds_simplify` (eval combines same-base powers + additive cancellation to expose u-freeness).
2. `nth_algebraic` recurses a Lagrange equation's radical branches (`y=2xy'+y'²` →
   `y'=−x±√(x²+y)`) through every early method → polyshift ran on them and hung the suite;
   fixed by lesson 1. Test early methods on nth_algebraic/Factorable branches, not just targets.

All DSolve ctest suites (dsolve_tests + m5/m12/m14/m17/m18/m19/m20 stress) + check-c99 green.
Corpus delta + STATUS/gate-baseline pending the full re-run (0-regression check across all
first-order solves, since early-polyshift placement affects them all).
