# M27 — §2.2.7 corpus (Problems 601–700) + corpus-harness system verification

Plan: `/Users/user/.claude/plans/let-s-continue-our-implementation-lively-ladybug.md`

## Phase A — §2.2.7 corpus file
- [x] Fetch `indexsubsection16.htm` + convert → `DSolve_test_status/DE_examples_227.m` (100 recs)
- [x] Verify §2.2.1–§2.2.6 regenerate byte-for-byte identical records (converter change record-neutral)
- [~] Subscripted-system indvar polish — SKIPPED (functional as-is; #609 solves; not needed)

## Phase B — Harness system verification (CORE)
- [x] Generalized `dsExplicitQ` (dsRuleFn/dsFnList) to accept List `fn`
- [x] Un-skip systems in `dsolveCheckCode`; route through DSolve + branch verifier
- [x] Header docstring updated (fn-list exclusion unneeded: residual free of fns post-subst)
- [x] Runner needs no change (verified by running §2.2.7/§2.2.6)

## Phase C — Report bucketer + re-baseline prior sections
- [x] `tools/dsolve_corpus_report.py`: scalar + systems lines in overall report
- [x] Re-baseline §2.2.6 (26 systems) → 95/100, 0 FAIL, baseline 2→5; reports + CMake + STATUS
- [x] Re-baseline §2.1.2 (204 systems) — DONE: 553/1204, 0 FAIL. Scalar 446/1000 (+10), systems 107/204. Baseline 572→655; reports/2.1.2.{tsv,md} + STATUS + CMake updated.
- [x] Triage FAILs — §2.2.6: 0 new FAIL; §2.2.7: fixed #684 (Solve collision) + #695 regression

## Phase D — §2.2.7 wire-up + chase coverage
- [x] Added `dsolve_corpus_2_2_7_tests` to `tests/CMakeLists.txt` (baseline 7)
- [x] Baseline-measured → reports/2.2.7.{tsv,md}. 90→93/100 after 2 engine fixes.
- [x] Chased misses: fixed separated-exp Simplify hang (+2 sys), Solve collision (+1). #604/#608 documented.

## Phase E — Docs + anti-overfit units
- [x] STATUS.md §2.2.7 block + M26/M27 wave-history + §2.2.6 re-baseline (§2.1.2 pending numbers)
- [x] README.md Contents row + systems-now-verified note
- [x] DSOLVE_PLAN.md M27 entry
- [x] docs/spec/changelog/2026-09-07.md summary + POSSIBLE_ZEROQ_IMPROVEMENTS.md #2
- [x] `t_m27_*` units in tests/test_dsolve.c (all pass)

## Verification
- [x] make + make check-c99 green
- [x] dsolve/solve/reduce/solve_corpus + all dsolve stress suites green (no regression)
- [x] dsolve_tests (incl t_m27_*) green
- [x] §2.1.2 (655) + §2.2.6 (5) re-baselined; §2.2.7 (7) 0-FAIL; §2.2.1–5 gates hold (2.2.1/2.2.2 improved)
- [x] valgrind spot-check: no new leaks (delta = documented per-call Integrate/Solve leak)
- [~] code-review-graph rebuild — MCP server was DOWN at session start (CONNECT_TIMEOUT); could not

## Review

**M27 complete.** §2.2.7 (Problems 601–700) added to the corpus, and the corpus harness
now VERIFIES systems (the headline change) rather than skipping them.

Results (all 0 FAIL):
- §2.2.7: **93/100** — 48/50 scalar + 45/50 systems (baseline 7).
- §2.2.6 re-baseline: **95/100** — 72/74 scalar + 23/26 systems (baseline 2→5).
- §2.1.2 re-baseline: **553/1204** — scalar 446/1000 (+10 from the Solve fix), systems
  107/204 (baseline 572→655).
- §2.2.1–§2.2.5 gates hold; §2.2.1 (4→2) and §2.2.2 (8→7) improved as a side-effect.

Three fixes (all root-cause, all verified):
1. Harness system verification (`dsolve_corpus_prelude.m`): `dsExplicitQ` on a List fn.
2. Separated-exponent `Simplify` hang (`dsolve_linsys.c`): tidy Expands exponential bodies.
3. Solve periodicity-index collision (`solveinv.c` + `dsolve_common.c`): fresh mint index +
   Element[…,Integers]-scoped family collapse. Fixed the M21-flagged wrong-answer class.

No regressions: dsolve/solve/reduce/solve_corpus + 7 dsolve stress suites + check-c99 green;
valgrind leak-clean. 3 anti-overfit units `t_m27_*`. Docs: DSOLVE_PLAN M27, STATUS, README,
changelog, PZQ backlog #2.
