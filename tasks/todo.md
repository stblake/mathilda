# Minimize campaign II — gaps from the 21–40 stress suite

Plan: `~/.claude/plans/here-are-20-more-flickering-tide.md`. Scope: everything
achievable (6 phases). Each phase: implement → test → re-run 21–40 → version
bump + tag → changelog + calculus.md → `make check-messages`/`check-c99`.

Measured baseline (v0.312): 4 solved (22,27,35,38), 4 timeout (24,31,32,34),
12 decline. Zero wrong answers.

## Phases
- [x] **P1 — Separable/additive decomposition (unconstrained)** → #24 SOLVED
      (60s→16ms). `mz_separable_min` + `mz_unconstrained_min`. minimize_tests green.
- [x] **P2 — General compact-region certificate skip** → #23 SOLVED (45ms), #34
      SOLVED (0.43s, -8/243). Generalized Fritz-John singular enumeration
      (`mz_add_singular_candidates`), single-probe `Resolve[Exists]` bounded check,
      shortcut-before-certificate. minimize_tests green.
- [x] **P3 — Equality-constraint variable elimination** → #31 SOLVED (60s→2ms).
      `mz_eliminate_solve`: bare-var equality substitution + reconstruction.
      minimize_tests green.
- [x] **P4 — Rational-function objectives** → #21 SOLVED (24ms), trivial rationals
      SOLVED; #40 soundly declines (boundary poles). `mz_rational_solve`:
      w=p/q auxiliary-variable reformulation + strict-sign gate. minimize_tests green.
- [x] **P5 — Integer: parametric Diophantine + coercive-integer** → #37 SOLVED
      (0.53s, {3,{-1,1,1}}). `mz_integer_parametric` (family→ellipsoid box) +
      `mz_strip_domain_element` (full-coverage Element→Integers; #25 mixed stays
      declined, sound). minimize_tests green.
- [x] **P6 — Constrained positive-dimensional QE (scaffold)** → settles
      `Minimize[{x+y,x<=y^2}]` = -Infinity (was decline). `mz_qe_infimum(cons)` +
      `mz_qe_witness(cons)` with attainment gate. #28/#36 stay declined (CAD/QE
      engine can't eliminate yet); #33 declines fast at the non-poly gate.

## Deferred (sound declines, out of scope): #25 mixed-integer, #29/#30/#39
transcendental, #32 multivariate radical sum.

## Review (v0.318, campaign complete)

**Scorecard 21-40: 4→10 solved, 4→2 timeout, 12→8 decline, zero wrong answers.**
- New solves: #21 (rational), #23 (Viviani compact), #24 (separable), #31
  (elimination), #34 (compact simplex sextic), #37 (integer parametric). Plus
  `{x+y,x<=y^2}` → -Infinity.
- Remaining declines all sound: #25 mixed-int, #28/#36 (CAD/QE-blocked), #29/#30
  transcendental, #33 multivariate Max/Abs, #39 transcendental const E, #40
  rational boundary poles.
- Remaining timeouts: #32 (multivariate radical sum, unchanged from baseline),
  #26 (6-var paraboloid↔plane).
- **Regression noted:** #26 went from a 9ms fast-decline to a ~60s slow-decline
  (still sound/terminates). Cause: Phase 3 elimination (z→x²+y²) exposes a 5-var
  residual whose internal Solve/Reduce probe overruns TimeConstraint — the
  pre-existing "TimeConstrained can't interrupt tight C loops" limitation. Not a
  correctness issue; flagged for the separate TimeConstrained-interruption fix.

**Verification:** `tests/test_minimize.c` +8 campaign tests (each win + a decline
soundness pin), all green. Build clean (gcc -O3 -Wall -Wextra, no warnings).
`make check-messages` ✓, `make check-c99` ✓. valgrind: **0 leaks allocated in
minimize.c**; ~20 extra blocks are the pre-existing rat_canon/CAD leak
(tasks/flint_ratcanon_leak.md) exercised by the extra Reduce calls.

**Docs:** `docs/spec/changelog/2026-10-05.md` (v0.313–v0.318 section),
`docs/spec/builtins/calculus.md` (new capabilities + trimmed Deferred),
`src/version.h` → 0.318.

**Pending user go-ahead:** git commit + per-phase tags (v0.313–v0.318).
