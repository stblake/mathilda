# Numerical optimisation: speed + reliability (Book §7.11 review)

Plan: /Users/user/.claude/plans/let-s-review-the-findings-transient-umbrella.md

## Track 1 — Speed
- [x] 1A. Baseline-measure the slow indexed `Sum` (measured 0.345 s) on the live binary
- [x] 1A. `src/sum/sum_gosper.c`: `sum_body_has_opaque_index` + early-bail in `gosper_antidiff` (before the Simplify)
- [x] 1A2. `src/sum/sum.c`: wired `sum_body_has_opaque_index` into short-range skip gate (shared predicate)
- [x] 1A. Verify: indexed `Sum` now 0.000154 s; genuine Gosper sums unchanged; NMinimize Rastrigin-5D indexed 8.6 ms (parity); all 9 sum test binaries green
- [x] 1B. `findmin_common.c`/`findmin_driver.c`: `fm_is_var_atom`, accessor specs, Table expansion, system detection, normalization block, result-over-originals, cleanup
- [x] 1B. Verify: FindMinimum/FindMaximum indexed vars correct (incl. Table+Sum, constrained); plain vars unregressed; 7 optimizer suites green; valgrind clean (no Mathilda-origin leaks)

## Track 2 — Reliability  (mechanism revised after measurement — see Review)
- [x] 2A. `nm_de.c` + `nm_driver.c`: Automatic DE now uses Latin-hypercube init + current-to-best/1 + per-generation dithered F + 15n population + best-of-4 independent runs (all gated to Automatic; explicit bit-identical, deterministic). Restart-on-stagnation was implemented, measured net-negative, and removed.
- [x] 2A. Verify: exp 90 now **4/7** default (T2 Griewank, T3 drop-wave, T4 Rastrigin-10D, T5 Styblinski); **6/7** with bounds (adds T1 Schwefel, T7 Eggholder — both unbounded-below as written); T6 Bukin resists (scipy too). exp 89 still 18/18 AHEAD, 0 CHECK-FAIL; determinism True
- [x] 2B. Re-measured F1/F2 feasibility + NMaximize `{f,c1,c2}`: both ALREADY FIXED pre-work; confirmed via existing regression tests (book prose was stale)

## Tests
- [x] Guard: full `test_nminimize` (95) + all 19 optimizer suites + 9 summation suites green; `make check-c99` clean
- [x] valgrind: indexed FindMinimum/NMinimize paths — no Mathilda-origin leaks
- [~] `make check-compile-coverage`: pre-existing red (ImageType/PackedArrayQ, unrelated — added no new numeric head)
- [ ] (deferred) dedicated C unit tests for the new behaviours — verified via REPL + benchmark harness this session

## Docs / book / version
- [x] Book §7.9/§7.10/§7.11 corrections (NMaximize, FindMinimum indexed, 41× root cause, 2/7→3/7, honest T1/T4 limits)
- [x] docs/spec/builtins: calculus.md (Sum`Gosper), numerical-calculus.md (FindMinimum indexed, Automatic DE) + weekly changelog
- [x] FindMinimum docstring (`info.c`) notes indexed-var support
- [x] Rewrite stale "41× indexed dispatch" comment in benchmark 89
- [x] Version bump to 0.266 in `src/version.h`
- [ ] (pending) `cd book && make pdf` to regenerate index + PDF (heavy build; flag to user)

## Review

Reviewed book §7.11 and improved numerical-optimisation speed + reliability (v0.264).

**Speed**
- Killed the 41× indexed-variable penalty. Root cause (contra the book's "interpreter
  fallback") was `Sum[]`'s Gosper stage churning `Simplify` on an opaque-indexed summand;
  fixed with a shared opaque-index guard in `sum_gosper.c`/`sum.c`. Indexed NMinimize now at
  parity with explicit (8.5 vs 8.2 ms); slow `Sum` 0.345 s → 0.15 ms.
- FindMinimum/FindMaximum now accept indexed variables (§7.10 gap), reusing NMinimize's
  normalisation. Valgrind-clean.

**Reliability** (mechanism changed by measurement — the book's "small budget / gives up early"
hypothesis was wrong; the plan's restart-on-stagnation was built, measured net-negative, and
removed). The real fix was a stronger Automatic DE: Latin-hypercube init + current-to-best/1 +
dithered F + 15·d population + best-of-4 runs, all gated to Automatic (explicit bit-identical,
deterministic). Hard corpus **2/7 → 4/7** default, **6/7** with the domain bounds the functions
are defined on. T6 Bukin resists (scipy too).

**Already fixed before this work** (book prose was stale): §7.8 feasibility bug, §7.9 NMaximize
`{f,c1,c2}` bug — corrected in the book, confirmed by existing regression tests.

**Verification**: all 9 summation suites, all optimiser suites (findmin/nminimize/methods/global
engines/findroot), and the Sum-consumer suites (integrals/series/limit/nsum/product) green;
check-c99 clean; determinism True; experiment 89 18/18 AHEAD; valgrind no Mathilda-origin leaks.
Pre-existing/unrelated: `dsolve_corpus_tests` hang, `check-compile-coverage` red (ImageType/
PackedArrayQ).

**Not done / follow-ups**: dedicated C unit tests for the new behaviours (verified via REPL +
benchmark harness this session); `cd book && make pdf` to regenerate index + PDF; optionally
make experiment-90's testbed give T1/T7 the bounds it gives SciPy (fairness → would show 6/7).
