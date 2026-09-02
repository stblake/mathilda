# DSolve M10 — Lie `bivariate` heuristic

Scope this pass: implement **`bivariate`** only (degree-bounded polynomial symmetry),
the exact generalization of the working `lie_linear`. Then reassess.

## Plan — DONE (bivariate landed; reassess remaining heuristics next)
- [x] Add `lie_poly_build` (build Σ coeff·x^i·y^j, i+j≤d) + `lie_coeff_name` helpers
- [x] Add `lie_poly_symmetry(omega, degree, ...)` — generalized `lie_linear` (NullSpace)
- [x] Add `lie_bivariate` (tries d=2 then d=3); chain into `dsolve_lie_try`
- [x] Update `dsolve_lie.c` header staging comment + docstring parenthetical
- [x] Unit tests (BV1, BV2, alias) `t_lie_bivariate` in `tests/test_dsolve.c`
- [x] Stress family `t_stress_lie_bivariate` in `tests/test_dsolve_stress.c` (8 A(u) members)
- [x] Build clean (`-std=c99 -Wall -Wextra`); all 3 `ctest -R dsolve` green
- [x] `make check-c99` exit 0
- [~] valgrind per-call: bivariate's own code is leak-clean, BUT solving triggers a
      PRE-EXISTING FLINT-bridge/rat_canon leak (56 B/call) — filed
      `tasks/flint_ratcanon_leak.md`; deferred per user decision (ship + track)
- [x] REPL spot-check + `?DSolve`LieSymmetry`
- [x] Docs: DSOLVE_PLAN.md, design doc §4, builtins/calculus.md, changelog
- [ ] Rebuild code-review graph
- [x] Reassessed with user → SHIP bivariate; track the leak separately

## Reassess next pass (roadmap in the plan file)
Remaining Lie heuristics: `abaco1_product`, `function_sum`, `abaco2_similar` (L2);
`chi`, `abaco2_unique_unknown`, `abaco2_unique_general` (L3). Math + forward-generated
test ODEs for all six already validated (see plan file "Next increments").

## Test ODEs (forward-generated, CORRECTED — verified bivariate-only in REPL)
The Plan agent's E1/E2 (`y/x+y²/x³`, `x²y/(x³+y²)`) turned out to have AFFINE symmetries
(E1: ξ=x,η=2y; E2: ξ=2x,η=3y) → caught by `lie_linear`, so they pass VACUOUSLY. Rejected.
Corrected isolation ODEs (deg1 nullspace = 0, abaco1_simple A/B/C all False, deg2 ≥ 1;
verified they solve via bivariate with residual PossibleZeroQ=True):
- **BV1**: `y' = -1/x + y/x + y^2/x^3`  → `-1/x + (1/2)Log[x+y] - (1/2)Log[-x+y] == C[1]`
- **BV2**: `y' = y/x^2 + y/x + y^2/x^3`  → `-1/x - Log[y] + Log[x+y] == C[1]`
- (spare) `y' = y/x + y/(x(x+y))` → `1/x - Log[x] + y/x + Log[y] == C[1]`
Symmetry for the A(u)-family `ω = y/x + A(y/x)/x` is ξ=x², η=xy.
Regression guard: `y'=y^2+x` (Riccati, existing `t_lie_declines`) STILL declines with
bivariate on — confirmed (its deg-2 candidate is tangent / non-elementary).

## Notes
- Safety-net: `lie_check` (†) gate + `dsolve_verify_implicit` back-sub ⇒ a mis-derived
  ansatz declines, never wrong. Bivariate reuses `lie_linear`'s exact pipeline.
- Reference: `benchmarks/.venv/.../sympy/solvers/ode/lie_group.py:478` (lie_heuristic_bivariate).

## Review — abaco2_similar increment (2026-09-02, committed 574d92a4, pushed to main)
- [x] `abaco2_similar` (§4.3) landed: [F(x),H(x)]+inverse; first heuristic to reach
      IRRATIONAL ω — solves y'=Sqrt[a x+b y+c], (a x+b y+c)^p. Cascade: after
      abaco1_product, before bivariate.
- [x] Unit `t_lie_abaco2_similar` + stress `t_stress_lie_similar` (8 radical ODEs, isolating).
- [x] **NthAlgebraic Root-guard** (dsolve_nth_algebraic.c): declines implicit-Root[] branches
      → fixes production hang on `y==x(y')²+(y')³` (Lagrange now claims it); explicit-root
      cases still solve.
- [x] **init.m loaded in all 3 DSolve test mains** (test_load_init_m in test_utils.h) — the
      env divergence that hid the hang. Watchdog 60s→120s. See memory
      project_dsolve_tests_load_init_m.
- [x] **Simplify power-distribution leak fixed** (simp_power.c, 3 branches) — pre-existing,
      surfaced by the Lie ds_simplify finalize; expr_new_function copies the args array so the
      heap holder must be freed.
- [x] Verified: 3 DSolve suites + 10 simplify/power suites + check-c99 green; per-solve valgrind
      leak-flat (1× vs 6× distinct calls identical after the simp_power fix).
- Remaining Lie heuristics: `function_sum` (L2); `chi`, `abaco2_unique_unknown`,
  `abaco2_unique_general` (L3).
