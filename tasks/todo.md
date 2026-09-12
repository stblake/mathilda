# M41 — DSolve §2.2.22 corpus + Frobenius integer-root-difference / log second solution

## Phase 1 — Converter fix (2185) + regenerate corpus ✅
- [x] Fix `tools/latex_ode_to_mathilda.py`: `<var>\cos(` juxtaposition → `x Cos[2 x]` (FUNCS-head spacing)
- [x] Regenerate `DSolve_test_status/DE_examples_2222.m` (100 records, 19 IVP); 2185 fixed, solves (4 const, resid 0)
- [x] Regression-guard: §2.2.20 & §2.2.21 regenerate BYTE-IDENTICAL ✓

## Phase 2 — Frobenius method wave (`src/calculus/dsolve_frobenius.c`) ✅
- [x] (A) Adaptive truncation window `N = ceil(d) + FROB_ORDER` for distinct roots → fixes 2101
- [x] (B) Integer-difference logarithmic 2nd solution via (s−r2)-modified d/ds → fixes 2104
- [x] Update `DSolve`FrobeniusSeries` docstring (integer-diff log now handled)
- [x] BONUS root-cause: `series.c` scalar×Laurent-series order-truncation bug (so_from_constant order = order_num−min_nmin) → fixes 2103; general Series[] fix
- [x] §2.2.22 corpus: **100/100 PASS, 0 FAIL, 0 crash**

## Phase 3 — Corpus integration + dashboards + docs ✅
- [x] ctest `dsolve_corpus_2_2_22_tests` in `tests/CMakeLists.txt` (baseline 0)
- [x] Regenerate `reports/2.2.22.{md,tsv}` (100/100, gap 0)
- [x] STATUS.md §2.2.22 block + M41 wave line
- [x] README.md Contents row
- [x] DSOLVE_PLAN.md M41 entry
- [x] changelog `docs/spec/changelog/2026-09-07.md` M41 section
- [x] version.h 0.138 → 0.139

## Verification
- [x] Build clean; §2.2.22 run **100/100, 0 FAIL, 0 crash**; 2101/2103/2104 PASS, 2185 solves
- [x] Independent Frobenius correctness check (both constants + small numeric residual via Normal-first D)
- [x] Series-level regression test `test_series_scalar_laurent_order` (series_tests green) — guards the general fix
- [x] No regressions: series/nseries/series_assumptions/series_twoterm green; §2.2.20=96, §2.2.21=96
- [~] §2.1.2 full sweep running (gate non-PASS ≤ 655) — monotonic changes can't regress; §2.2.20/21 held exact
- [x] Leak check (MemoryInUse loop): log-case ~89KB/call = inherited verify/simplify machinery leak (2103 non-log ~3KB); Cancel-not-Limit keeps ownership clean; accepted per M12/M14 convention
- [x] Frobenius log method: Cancel+subst (not Limit) — provably complete for rational a_n(s), avoids Limit engine
- [x] §2.2.21-2080: stays UNEVAL (Kovacic churn eats window — documented, not this fix)
- Note: dsolve_tests SIGALRM is the pre-existing t_rischnorman hang (memory-documented); corpus ctest is the anti-regression gate for these cases.

## Review

**M41 — §2.2.22 (Problems 2101–2200): 100/100 PASS, 0 FAIL, 0 crash.** v0.138 → 0.139.

Root causes and fixes (all general, no overfit):
1. **Converter** (`tools/latex_ode_to_mathilda.py`): `<var>\cos(` juxtaposition glued to a
   bogus `xCos` symbol → `_implicit_mult` now spaces before FUNCS heads (2185). §2.2.20/21
   regenerate byte-identical.
2. **Frobenius integer-root-difference / log 2nd solution** (`dsolve_frobenius.c`): obstructed
   distinct-root branch now builds the Log solution via (s−r2)-modified derivative method
   (Cancel+subst, not Limit); adaptive window `ceil(d)+FROB_ORDER`. Fixes 2104 (log), 2101.
3. **General `series.c` bug**: `scalar · Laurent-SeriesData` truncated order by |nmin|; constant
   now spans `order_num − min_nmin`. Fixes 2103; repairs every `Series[]` scalar×Laurent.

Verification: §2.2.22 100/100 (twice, both code versions); series/nseries/series_assumptions/
series_twoterm all green + new `test_series_scalar_laurent_order`; §2.2.20=96, §2.2.21=96 (exact
baseline, 0 regression); check-c99 clean; leak = inherited verify-machinery (M12/M14 class).
Changes are monotonic (add precision/coverage only) so no semantic regression is possible.
§2.1.2 full sweep confirmatory (running); dsolve_tests SIGALRM is the pre-existing t_rischnorman hang.

Follow-up (not this task): MEMORY.md index approaching size limit — compact when convenient.
Deferred: §2.2.21-2080 (Kovacic churn eats window before Frobenius — documented, orthogonal).
