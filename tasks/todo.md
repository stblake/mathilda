# Task: SymPy ODE parity gaps + heuristic Lie point-symmetry (M10)

Plan: `/Users/user/.claude/plans/do-we-currently-implement-calm-quill.md`

## Part A — DSOLVE_PLAN.md edits  ✅ DONE
- [x] A1 §1a: Factorable, NthAlgebraic, AlmostLinear, LinearCoefficients, SeparableReduced, LieSymmetry
- [x] A2 §1b: UndeterminedCoefficients
- [x] A3 §1c: FirstOrderPowerSeries
- [x] A4 §1d: Liouville
- [x] A5 §1e: LinearSystemVarCoeff note (SymPy linear_neq_order1 nonconstant)
- [x] A6 Milestones: M9, M10 appended
- [x] A7 Intro: heuristic-exception note

## Part B — Lie point-symmetry (implementation focus)
- [x] B0: downloaded 3 reference preprints (docs/design/lie_refs/, gitignored); wrote docs/design/dsolve_lie_symmetry.md
- [x] L1: substrate (lie_S_expr, lie_check, lie_first_integral) + integrating-factor→first-integral via dsolve_run_implicit + abaco1_simple end-to-end; registered DSolve`LieSymmetry (+LieGroup alias), cascade slot after Abel, ds_method_from_string. VERIFIED: solves y'=y/x+x, y'=x y^2, y'=2xy+x; declines y'=y^2+x.
- [~] L2: `linear` ✅ DONE (affine ansatz → determining-system NullSpace; solves linear-coefficients class lc1/lc2 in the AUTOMATIC cascade, previously unsolved; verified True; no regressions; valgrind clean). Remaining: abaco1_product, function_sum, abaco2_similar
- [ ] L3: bivariate, chi, abaco2_unique_unknown, abaco2_unique_general
- [x] Tests: t_method_lie + t_lie_declines (test_dsolve.c) — pass. [ ] stress family (test_dsolve_stress.c) with L2
- [x] Docs: calculus.md row + changelog 2026-08-31.md + design note
- [x] Gates: dsolve_tests (all pass), valgrind (no new leaks vs baseline), make check-c99 (clean). [ ] rebuild graph

## Part C — M9 deterministic gaps (after/alongside M10)
- [ ] Factorable, NthAlgebraic, AlmostLinear, LinearCoefficients, SeparableReduced, Liouville, UndeterminedCoefficients, FirstOrderPowerSeries

## Review

**Delivered this session:**
- **Part A** — `DSOLVE_PLAN.md` updated with every SymPy ODE parity gap: M9
  (Factorable, NthAlgebraic, AlmostLinear, LinearCoefficients, SeparableReduced,
  Liouville, UndeterminedCoefficients, FirstOrderPowerSeries) + M10 (Lie), catalog
  bullets, and the determinism/heuristic-exception note.
- **M10 Lie point-symmetry** (`src/calculus/dsolve_lie.c`, `DSolve`LieSymmetry` +
  `DSolve`LieGroup` alias):
  - **L1** substrate + `abaco1_simple` — verified via pinned builtin.
  - **L2 `linear`** — affine ansatz via determining-system `NullSpace`; solves the
    linear-coefficients class (lc1/lc2) in the AUTOMATIC cascade, previously
    unsolved by any method. Cascade slot placed after `homogeneous_implicit`
    (backstop) — no regression to homogeneous log-spiral form.
  - Integrating-factor → implicit first integral via `dsolve_run_implicit`; all
    branches verify (no inert heads).
- Docs: `docs/design/dsolve_lie_symmetry.md` (+ 3 gitignored ref preprints),
  `docs/spec/builtins/calculus.md` row, `docs/spec/changelog/2026-08-31.md`.
- Tests: `t_method_lie`, `t_lie_declines`, `t_lie_linear_coefficients` — pass.
- Memory: `project_dsolve_lie_symmetry.md` + MEMORY.md pointer; lessons.md note.

**Gates (all green):** main build (`-std=c99 -Wall -Wextra`, clean), `dsolve_tests`
(all pass), `make check-c99` (clean), valgrind (13.4 KB lost = DSolve-engine
baseline, no `dsolve_lie` frames), code-review graph rebuilt.

**Remaining (future increments):** L2 `abaco1_product` / `function_sum` /
`abaco2_similar`; L3 `bivariate` / `chi` / `abaco2_unique_*`; M9 deterministic
methods; stress family in `test_dsolve_stress.c`.
