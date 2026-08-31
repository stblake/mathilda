# DSolve M9 — Test Hardening (tests + docs only; no src changes)

Make DSolve's test suite extensive and uniform: backfill thin unit coverage,
add pinned-method tests where a builtin exists, and add forward-generator stress
families for every method group that lacks one. Verify every case against the
current binary before committing (no vacuous or declining tests).

## Phase 0 — Validate case corpus against ./Mathilda
- [x] Confirmed backfill unit cases (Homogeneous, ReductionOfOrder, Clairaut, PowerSeries)
- [x] Confirmed every stress-family generator stays in-domain across its grid

## Phase A — Backfill unit tests in tests/test_dsolve.c
- [x] Homogeneous: t_method_homogeneous + t_homogeneous_more + t_ivp_homogeneous
- [x] ReductionOfOrder: t_method_reduce_order + t_reduce_order_more
- [x] Clairaut: t_method_clairaut + t_clairaut_more
- [x] PowerSeries: t_powerseries_more (2 ordinary forms)
- [x] Registered all 8 new t_* in main()

## Phase B — New stress file tests/test_dsolve_stress.c
- [x] Forward generators: LinearFirstOrder (target-solution), Separable, Bernoulli,
      Homogeneous (curated), Exact (from potential), ConstCoeff (from spectrum),
      Euler (from indicial), ReductionOfOrder
- [x] Systems 2x2 (distinct/complex/defective/singular) + triangular varcoeff
- [x] PDE first-order (transport/forcing) with C[1][z_]:>Sin[z]
- [x] Runs in ~3s (well under 60s alarm) — no split needed

## Phase C — Wiring
- [x] add_executable/add_test for dsolve_stress_tests in tests/CMakeLists.txt

## Phase D — Build, run, gates
- [x] cmake + make dsolve_tests dsolve_m5_stress_tests dsolve_stress_tests
- [x] ctest -R dsolve --output-on-failure (3/3 green)
- [x] make check-c99 (exit 0)
- [x] valgrind: definitely lost == 13,608 bytes == engine baseline (not per-call)
- [x] rebuild code-review graph

## Phase E — Docs
- [x] docs/spec/changelog/2026-08-31.md: "DSolve — test hardening (M9)" section
- [x] DSOLVE_PLAN.md Testing section updated
- [x] tasks/lessons.md: target-solution generator + engine-baseline lessons

## Review

Delivered tests-only DSolve hardening (no `src/` changes):

- **Part A** — 8 new unit tests in `tests/test_dsolve.c`: pinned-method coverage
  for the three scalar methods that lacked it (Homogeneous, Clairaut,
  ReductionOfOrder) + more in-domain forms and an IVP for the four thin methods.
- **Part B** — new `tests/test_dsolve_stress.c` (`dsolve_stress_tests`, ~75
  generated cases, ~3s): 11 forward-generator families spanning the rest of the
  cascade. Every case guards `Head === List` before the residual check.
- Gates: 3/3 dsolve ctest targets green; `check-c99` clean; valgrind leak equals
  the pre-existing 13.6KB one-time engine baseline.

**Solver gaps surfaced (not fixed — out of scope, reported):**
1. `DSolve`Homogeneous` needs the separated integral to invert explicitly and
   declines otherwise (e.g. `y'==(x+y)/(x-y)`, `y'==(x+2y)/(2x+y)`) — only a
   subset of rational-homogeneous RHS solve, hence the curated Homogeneous family.
2. `DSolve`ReductionOfOrder` declines `y''==1+(y')^2` and `y''==-2x(y')^2` (first-
   order-in-p sub-solve / inversion) — the `f(x)·y'` family stays in-domain.
3. `y'+x·y==Sin[x]` declines (integrating factor `Exp[x²/2]`, non-elementary
   `∫Sin·Exp[x²/2]`) — a genuine elementary-integrability boundary, handled by the
   target-solution generator rather than a fix.

---

# Follow-up: fix Homogeneous (gap 1) + ReductionOfOrder (gap 2)  ✅ DONE

User asked to address gaps 1 and 2 above. Root-caused each, then a minimal
`src/` fix per gap (correctness gated by the existing back-substitution verify).

- **Gap 2 — `dsolve_reduce_order.c`.** Root cause: the `D[∫p]==p` guard used
  `zero_test_decide`, too weak for a correct-but-unsimplified antiderivative
  (multi-`Log` `∫Tan`, `1/(C(1+x²/C))` from `∫1/(C+x²)`). Fix: escalate the guard
  to `PossibleZeroQ` (matches `integrate_derivdivides.c`); a wrong `Integrate->0`
  still fails as `-p`, and the final verify is the backstop. Opens `1+(y')^2`,
  `-2x(y')^2`, the autonomous `a+b(y')^2` and Riccati-in-p `c x (y')^2` families.
- **Gap 1 — `dsolve_homogeneous.c`.** Root cause: `Solve` cannot invert a
  sum-of-logs relation. Fix: new `homog_exp_log_invert` fallback — exponentiate to
  `Prod g_i^{c_i} == C[1] x`, clear fractional exponents by raising to power
  `d ∈ {1,2,3,4,6}` (`PowerExpand`), `Solve` → `Root` branches; `dsolve_run`
  filters spurious sheets. Opens the real-root rational-homogeneous family
  (`(x+2y)/(2x+y)`, …); the transcendental ArcTan subset still declines (no
  explicit inverse — would need implicit-solution support).

Tests: `t_homogeneous_algebraic`, `t_reduce_order_riccati` (unit) + broadened
Homogeneous stress + new `t_stress_reduce_order_riccati`. Gates: 3/3 dsolve ctest
green, no regressions; `check-c99` clean; valgrind shows **no per-call leak**
(1×==6× the algebraic solve == 13,440 bytes, the one-time engine baseline).

Remaining (reported, not done): the transcendental-implicit Homogeneous subset
(`(x+y)/(x-y)`) needs implicit-solution output — a larger substrate change.

---

# Follow-up 2: implicit first-integral solutions (transcendental homogeneous)  ✅ DONE

User asked to take on the implicit-solution subset. Added a parallel substrate
path rather than overloading the explicit `y[x] -> body` pipeline.

- **Substrate** (`dsolve_common.{c,h}`): `dsolve_run_implicit` +
  `dsolve_method_builtin_implicit`; a try-fn returns the LHS `G` of `G == C[1]`.
  `dsolve_verify_implicit` verifies by implicit differentiation (from
  `d/dx[G==C]=0`, `y' == -G_x/G_y` must satisfy the residual — same permissive
  not-decidably-nonzero policy as `dsolve_verify_body`). `dsolve_implicit_rhs`
  fits an IC as `G(x0,y0)`. Output `{{ G(x,y[x]) == C[1] }}` (an `Equal` branch).
- **`dsolve_homogeneous_implicit_try`** (`dsolve_homogeneous.c`): `G =
  (∫1/(F(v)-v) dv /. v->y[x]/x) - Log[x]`. Runs after explicit Homogeneous (both
  auto cascade — before series — and as the pinned-method fallback).
- Tests verify by implicit differentiation (the `Equal` branch can't use `/.`):
  `t_homogeneous_implicit` (unit, incl. IVP + pinned) + `t_stress_homogeneous_implicit`
  (5 transcendental forms). Updated the now-obsolete decline assertion in
  `t_homogeneous_algebraic` (the transcendental case solves implicitly now).

Gates: 3/3 dsolve ctest green, no regressions; `check-c99` clean; valgrind no
per-call leak on the implicit path (1×==6× == 13,440 + 6,312 bytes baseline).
Closes the last honest Homogeneous gap.
