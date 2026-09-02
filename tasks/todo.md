# Task: DSolve M9 — SymPy deterministic parity gaps (all 8 methods)

Plan: `/Users/user/.claude/plans/let-s-continue-our-dsolve-cheeky-music.md`

## Phases

- [ ] **P0** Baseline: confirm main build + `ctest -R dsolve` green before changes.
- [x] **P1** `Factorable` + `NthAlgebraic` (front; recurse+union) — DONE. New files
      `dsolve_factorable.c`, `dsolve_nth_algebraic.c`; cascade front slots; unit +
      stress tests; all dsolve ctest green. Key fixes: factor over plain-symbol
      substitutes with a `PolynomialQ` gate (FactorList HANGS/misfactors on raw
      funcapps); keep only DIFFERENTIAL factors (contain a derivative d_k, k>=1) so
      the trivial `y==0` branch never hijacks recursive callers (AutonomousReduction).
- [x] **P2** `LinearCoefficients` (`dsolve_lincoeff.c`, explicit + implicit two entries) —
      DONE. det!=0 shift->homogeneous (explicit Root-form / implicit log-spiral), det==0
      parallel v=a1x+b1y->separable (implicit). Superseded the Lie `linear` test
      (t_lie_linear_coefficients now pins DSolve`LieSymmetry). Unit + stress; suites green.
- [x] **P3** `AlmostLinear` (`dsolve_almostlinear.c`, u=Integrate[g,y]->linear, explicit) +
      `SeparableReduced` (`dsolve_sepreduced.c`, w=x^n y->separable, implicit) — DONE.
      Unit + stress; all suites green (no SR pre-emption regression).
- [x] **P4** `Liouville` (`dsolve_liouville.c`, 2nd-order two-quadrature) — DONE.
      g(y)(y')^2+h(x)y' form; EG==C[1]EH+C[2] solved for y. After autonomous in cascade.
      Unit + stress; suites green. (Cosmetic Solve::ifun on Log-inversion, MMA-consistent.)
- [x] **P5** `UndeterminedCoefficients` (`dsolve_undetcoeff.c`) — DONE. Exposed
      `dsolve_homog_basis` from constcoeff into common (shared, zero-regression). UC
      forcing (poly/exp/sin/cos + products + sums) via superposition + "increment s"
      resonance; before constcoeff in cascade. Key fix: Expand (not Simplify) the
      residual before Coefficient[.,Cos[b x]] (Simplify rewrites Cos[2x]). Unit + stress.
- [x] **P6** `FirstOrderPowerSeries` (extend `dsolve_frobenius.c`) — DONE. Iterative
      Taylor coeffs about x0=0; PINNED-ONLY (not auto — a first-order ODE with no closed
      form stays unevaluated, matching MMA/SymPy opt-in). Muted the pole probe. Unit + stress.
- [x] **P7** Cascade wiring finalized; docs (`calculus.md` 8 rows) + changelog
      (`2026-08-31.md`) + `DSOLVE_PLAN.md` (M9 ✅, catalog `[✓]`); graph rebuild.

## Per-phase checklist (each method)
new `.c` (3-fn contract) → COMMON_SRC → cascade + enum + ds_method_from_string +
extern + init + pinned switch → unit tests (`t_method_*`, variant, `*_declines`,
auto) → stress family → build + `ctest -R dsolve` green → valgrind spot-check.

## Gates (before done)
- [ ] main build `-std=c99 -Wall -Wextra` clean
- [ ] `ctest -R dsolve` (dsolve_tests, dsolve_stress_tests, dsolve_m5_stress_tests) green
- [ ] `make check-c99` exit 0
- [ ] valgrind: no new leak beyond ~13.6 KB engine baseline
- [ ] docs/spec/builtins/calculus.md rows + changelog 2026-08-31.md notes
- [ ] code-review graph rebuilt; memory/lessons updated

## Review

**Delivered:** all 8 M9 methods — `Factorable`, `NthAlgebraic`, `AlmostLinear`,
`LinearCoefficients`, `SeparableReduced`, `Liouville`, `UndeterminedCoefficients`,
`FirstOrderPowerSeries`. New files `dsolve_{factorable,nth_algebraic,almostlinear,
lincoeff,sepreduced,liouville,undetcoeff}.c`; `FirstOrderPowerSeries` extends
`dsolve_frobenius.c`. Shared `dsolve_homog_basis` + `dsolve_extract_applied_bodies`
added to `dsolve_common.{c,h}`; `dsolve_constcoeff.c` switched to the shared helper.
Cascade + enum + `ds_method_from_string` + pinned switch wired in `dsolve.c`;
`COMMON_SRC` updated. Unit tests (2–3 per method + declines) and forward-generator
stress families added.

**Gates (all green):** main build `-std=c99 -Wall -Wextra` clean; `ctest -R dsolve`
(dsolve_tests, dsolve_stress_tests, dsolve_m5_stress_tests) 100%; `make check-c99`
exit 0; valgrind — 7/8 new methods per-call flat at the 13,440 B one-time baseline.

**Known issue (pre-existing, out of scope):** `AlmostLinear` inherits a per-call
leak in the `Integrate` engine (`integrate.c` ~920/256) that `LinearFirstOrder` and a
bare `Integrate[Exp[x] Sin[x], x]` show identically (≈168 B/call). Flagged for a
separate core-integrator fix; NOT introduced by M9.

**Lessons captured:** FactorList hangs/misfactors on raw funcapps (factor over plain
symbols + PolynomialQ); differential-factor gate for recursion safety; Expand (not
Simplify) before Coefficient[·,Cos[b x]]; the Integrate leak. Memory + MEMORY.md
updated.
