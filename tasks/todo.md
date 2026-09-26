# Task: Stress-Testing `Refine`

Plan: `/Users/user/.claude/plans/let-s-thoroughly-stress-test-sorted-rain.md`

## Phase A — Corpus (159 cases) ✅
- [x] Verify harness mechanics (`./Mathilda -file`, timeout, SameQ probe)
- [x] Build `corpus.py` — 12 categories, first-principles expected
- [x] Sanity-check corpus (calibration probe surfaced gaps + wins)

## Phase B — Harness ✅
- [x] `run_stress.py` — per-case isolated process + hard timeout
- [x] Verdict logic: PASS / UNCHANGED / DIVERGENT / HANG / CRASH via SameQ
- [x] Full run → results JSON

## Phase C — Gap report ✅
- [x] `tasks/refine_stress_report.md` — matrix + prioritized failures + code locations

## Phase D — Safe fixes ✅ (v0.203)
- [x] G1 soundness: SZ sampler downgrades FALSE→UNKNOWN under coupling equality (zero_test.c) + Refine Equal TRUE via substitution (refine.c)
- [x] G3 Element Positive/Negative/NonNegative/NonPositive domains (simp_builtins.c)
- [x] G2 NonPositive rewrites, G4 Arg under sign, G6 Abs[complex], G5 symbol-set restriction, Log identities (simp_assume_rewrite.c)
- [x] Regression tests: test_stress_fixes (24 assertions) in tests/test_refine.c
- [x] Re-run corpus 137→151 PASS, 0 DIVERGENT, zero regressions in refine/simplify/pzq/element/assuming
- [x] version.h 0.202→0.203 + changelog + spec updated
- [x] Deep gaps catalogued as follow-up in the report

## Verification ✅
- [x] `make -j` clean; `make check-c99` clean
- [x] CMake test build + all shared-engine suites green
- [x] valgrind: no leak attributable to changed files (only macOS objc/dyld baseline noise)
- [ ] git tag v0.203 (pending user approval to commit)

## Review
- Built a 159-case first-principles corpus (SameQ-based, printer-independent, per-case process
  isolation). Baseline 137 PASS. Found and fixed one **soundness bug** (assumption-aware zero test
  returned a wrong `False` under coupling equalities — poisoned both `Refine` and `PossibleZeroQ`)
  and six missing-feature families; final **151/159 PASS, 0 wrong answers**. 8 open follow-ups all
  return the input unchanged (never wrong): integer-linear trig, `Element[√2,Algebraics]`,
  deep-positivity for non-strict/compound-real bases, the >6-var CAD ceiling, and `Sign` of an exact
  algebraic zero. The methodology also caught two of my own wrong expected values and a harness
  accumulator/data-symbol collision.
