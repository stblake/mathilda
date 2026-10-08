# Task: Implement `Minimize` / `Maximize`

Plan: `/Users/user/.claude/plans/pasted-content-id-2dbb-we-should-golden-hopper.md`
Scope: Tier 1 (exact univariate poly) + Tier 2 (multivariate unconstrained + constrained poly over Reals).
Discipline: exact oracle backs every result; any "don't know" ⇒ decline (NULL). Inexact input ⇒ NMinimize.

## Phase 1 — Infra skeleton (build green)
- [ ] Intern `SYM_Minimize` / `SYM_Maximize` in `src/sym_names.h` + `src/sym_names.c` (3 edits each)
- [ ] Prototype `builtin_minimize` / `builtin_maximize` in `src/numerical_calculus/findmin.h`
- [ ] New `src/numerical_calculus/minimize.c` — stub returning NULL
- [ ] Register in `findmin_init()` with `ATTR_HOLDALL | ATTR_PROTECTED`
- [ ] Docstrings in `src/info.c` (next to NMinimize)
- [ ] Option tables in `src/options_builtin.c` (WorkingPrecision)
- [ ] `tests/test_minimize.c` stub mirroring `test_nminimize.c`
- [ ] `tests/CMakeLists.txt`: add minimize.c to COMMON_SRC + add_executable(minimize_tests)
- [ ] `make -j` clean; minimize_tests builds

## Phase 2 — Arg normalization + numeric fallback
- [ ] `mz_normalize_args` (f / {f,cons}, vars, dom, options) reusing fm/nm parsers
- [ ] inexact-input detection
- [ ] `mz_numeric_fallback` → synthesize+eval `NMinimize[...]`
- [ ] tests: inexact → numeric; unsupported exact → unevaluated

## Phase 3 — Tier 1 (a) univariate polynomial
- [ ] `mz_poly_tail` (degree parity + leading sign; bounded/unbounded/attainment)
- [ ] `mz_univar_poly`: D → solvepoly → `mz_candidate_min` (rru_sign_compare)
- [ ] `mz_result_exact / _unbounded` builders
- [ ] univariate headline tests pass

## Phase 4 — Maximize
- [ ] `builtin_maximize` via `-f` negation + value negation + tag swap
- [ ] max tests pass

## Phase 5 — Tier 2 (b) multivariate unconstrained
- [ ] gradient system via `Solve[{grad==0}, vars, Reals]` / reduce_zerodim; decline if positive-dimensional
- [ ] `mz_entails(A,P,dom)` = eval `Reduce[A && !P]` == False
- [ ] global lower-bound certificate
- [ ] tests (isolated min passes; positive-dimensional declines)

## Phase 6 — Tier 2 (c) constrained polynomial
- [ ] KKT/active-set enumeration; closure of region
- [ ] feasibility filter + candidate min + lower-bound certificate
- [ ] `_infeasible` / `_notattained` outcomes
- [ ] constrained headline tests pass

## Phase 7 — Docs, release, verification
- [ ] docs/spec/builtins/ (optimization doc) + docs/spec/changelog/2026-10-05.md
- [ ] version bump 0.303 → 0.304 (src/version.h), commit note, tag v0.304
- [ ] full test sweep + valgrind + check-messages + check-c99 + rebuild graph

## Review

**Status: complete (v0.304). All 26 `minimize_tests` pass; `nminimize_tests` /
`findmin_tests` still pass (no regression). Not yet committed (on `main`).**

Delivered:
- `Minimize` / `Maximize` builtins (`src/calculus/minimize.c` — a symbolic
  *calculus* module, registered by `minimize_init()` from `core_init`),
  Protected (non-Hold, faithful to Mathematica — deviated from the plan's
  `HoldAll`, see note below). Symbols, docstrings, option tables, registration,
  CMake, docs, changelog, version bump all wired.
- Tier 1 (a): exact univariate polynomial + tail-theorem bounded/unbounded/
  attainment.
- Tier 2 (b)/(c): multivariate unconstrained (gradient + Reduce certificate)
  and constrained polynomial over Reals (KKT/active-set + closure + certificate
  + attainment), unified in `mz_exact_poly`; `Inequality` chains expanded; LPs
  fall out as the degenerate case.
- Numeric fallback to NMinimize for inexact input; exact-undecidable declines.
- Infeasible -> `{Infinity, ...}`+`infeas`; unbounded -> `{-Infinity, ...}`+`natt`.

Verified against every supported Mathematica example in the request (exact
match), plus the deferred cases decline cleanly. check-messages / check-c99
pass. Valgrind: `minimize.c` adds **zero** leaked bytes (420-block
definitely-lost total is pre-existing one-time Solve/Reduce/FLINT init leakage,
identical for `Print[1]` and for 55 Minimize calls).

Design note / deviation from approved plan: the plan specified `ATTR_HOLDALL`;
I used `ATTR_PROTECTED` only. Mathematica's real `Attributes[Minimize]` is
`{Protected, ReadProtected}` (non-Hold), and the sibling `NMinimize` is
deliberately non-Hold. Non-Hold is more faithful AND simpler: the exact engine
receives already-evaluated symbolic args and works by substitution, needing no
Block-binding. Documented in `findmin.c` and the changelog.

Deferred (sound declines, next tiers): transcendental closed forms, parametric
`Piecewise`, positive-dimensional minimizer sets, general unbounded/natt via QE,
exact `Integers`/ILP, `MinValue`/`ArgMin`.
