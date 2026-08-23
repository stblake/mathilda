# Reduce — Extensive Stress Testing

Plan: `/Users/user/.claude/plans/per-the-dev-work-clever-moth.md`

## Layer 1 — Expanded C unit tests (`tests/test_reduce.c`)  ✅
- [x] Capture exact FullForm ground-truth for every new case from `./Mathilda` (112 cases)
- [x] all black-box groups (decides, equations + decline, real inequalities,
      linear systems, parametric systems + decline, integer/rational, decline
      soundness net, known limitations) — 127 assertions
- [x] white-box `test_wb_*` (atom emit, logic expand, unsupported) — 23 assertions
- [x] Build + run `reduce_tests` green (150 assertions, exit 0)

## Layer 2 — Semantic corpus verifier (new)  ✅
- [x] `tests/reduce_check_prelude.m` (sample-grid equivalence + witness oracle)
- [x] `tests/reduce_corpus.m` (66 records, all PASS)
- [x] `tests/test_reduce_corpus.c` (fork-per-case runner, mirror test_solve_corpus.c)
- [x] `tests/CMakeLists.txt` (add reduce_corpus_tests target + add_test)
- [x] Build + run `reduce_corpus_tests` → 0/66 non-PASS
- [x] Adversarial check: oracle rejects wrong verdicts, incompleteness, missing
      roots, strict/non-strict boundary errors, spurious roots (non-vacuous)

## Verification  ✅
- [x] valgrind on reduce_tests: no NEW leaks from test code or reduce engine
      logic (delta vs original = ~8 blocks in the pre-existing Together/rat_canon
      path + init-time interned symbols; ~420 blocks are macOS system-lib noise)
- [x] `make check-c99` clean (exit 0; changes are tests/ only, src/ untouched)
- [x] docs/spec changelog updated (2026-08-17.md)
- [ ] Rebuild code-review-graph

## Review
- **Layer 1**: `tests/test_reduce.c` 264 -> ~450 lines, 15 test functions, 150
  assertions. Every string pinned from the binary; soundness-verified by hand.
- **Layer 2**: semantic corpus (66 cases) verifies logical-formula outputs by
  sampling (grid equivalence + witnesses) rather than string match — robust to
  ordering/spelling/non-minimal forms, and proven to detect real soundness
  defects. Registered in CMake with `add_test`; ships green at baseline 0.
- **No engine code changed** — per "capture & flag": quirks are pinned + flagged.

## Phase-8 quirk fixes (requested 2026-08-23) — ✅ DONE
- [x] Q1: `Reduce[a x == 0, x]` → `(a!=0 && x==0) || a==0` (reduce_eq.c: is_zero
      base case for the vanished lower-order remainder)
- [x] Q2: Reals linear-equation systems back-substitute (reduce_fm.c: generalized
      lo==hi equation collapse + forward constant back-substitution) →
      `x+y==1 && x-y==3, {x,y}, Reals` gives `x==2 && y==-1`
- [x] Q3: parametric conditions print minimally (reduce_sys.c: sys_norm_condition
      + collect_leaf_params) → `1-a!=0`→`a!=1`, `-1+2a==0`→`a==1/2`
- [x] Regression tests: unit pins moved out of test_known_limitations (deleted)
      into test_equations / test_linear_systems / test_parametric_systems;
      corpus gained eq-param-ax0 + fm-eq-determined, dropped stale dec-ax0-gap
- [x] reduce_tests (149) + reduce_corpus_tests (67) green; solve_tests +
      solve_corpus (97) unregressed; valgrind clean (no new leaks); check-c99 0;
      changelog + solutions-of-equations.md updated

## Notes
- All three fixes are soundness-preserving (sampling corpus) and leave every
  previously-correct output unchanged (a x==b, quadratic split, triangular
  regions, x+y==1 && x>0 stays x>0 && y==1-x).
- Pre-existing (unrelated) leak in the `Together`/rat_canon path is surfaced by
  Reduce's atom canonicalisation; not introduced by any of this work.

## Findings (quirks discovered)
- `a x == 0` declines though `a x == b` solves (Mathematica: `x==0 || a==0`)
- Reals linear-equality systems return sound-but-non-minimal forms (no back-substitution)
- Parametric conditions print non-minimally (`1 - a != 0` vs `a != 1`)
