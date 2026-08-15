# L-BFGS-B for FindMinimum — implementation todo

Plan: `/Users/user/.claude/plans/let-s-compare-the-existing-harmonic-jellyfish.md`
Branch: `feat/lbfgsb-findminimum`

## Design deviations from plan (all verified, documented here)
- Bound handling: ACTIVE-SET projection (not the exact BLNZ generalized-Cauchy-
  point/subspace-min). Reaches the same optima; the bound-active discriminating
  case (Rosenbrock on x<=0.5 -> {0.25, x->0.5, y->0.25}) is EXACT. GCP left as a
  future refinement.
- Line search: robust strong-Wolfe with expansion + best-point fallback (tries
  alpha=1 first; the shared fm_line_search's 1/||d|| cap throttled L-BFGS).
- Convergence: projected-gradient inf-norm + machine-noise relative-f stall
  (removed the single-step displacement test, which pre-empted convergence on
  narrow valleys).
- Indexed vars: NOT added to FindMinimum. Large-n cases use programmatic scalar
  symbols Symbol["z"<>ToString[i]] (no driver surgery).
- Extended Rosenbrock is multimodal for n>=4 -> large-n scaling uses the
  UNIMODAL ill-conditioned quadratic instead (method-independent optimum).

## M1 core — DONE
- [x] Enum, parse strings, needs_grad, driver dispatch, fm_run_penalty case, MPFR nimpl fallback
- [x] fm_run_lbfgsb: two-loop + curvature-skip + active-set + robust Wolfe
- [x] All scalar/large-n/bound/general cases converge (verified)

## M2 (true GCP) — DEFERRED (active-set variant reaches same optima)

## Tests — DONE
- [x] tests/test_lbfgsb.c (26 tests, 7 groups) — ALL GREEN
- [x] Registered in tests/CMakeLists.txt
- [x] findmin/nminimize no regression

## Benchmark — DONE
- [x] benchmarks/64-lbfgsb-scaling/lbfgsb_scaling.{m,py} (10 cases)
- [x] Ran exp 64 vs scipy: CHECK-FAIL=0, INCOMPLETE=0; AHEAD 6 / SLOWER 4
- GOTCHA fixed: build objectives with Total[Table[...]] NOT Sum[...] — Sum
  attempts a closed form with a symbolic iterator (Symbol["z"<>ToString[i]]
  collapses to one concrete symbol -> geometric series) and burns the whole
  budget. Tests dodge it via v[[i]] (Part blocks the closed-form attempt).
- FINDINGS (honest, informative):
  * Within Mathilda: LBFGSB 402ms vs QuasiNewton 1650ms at n=1000 (4x — the
    scalability win); QN n=50->1000 grows 609x vs LBFGSB 128x.
  * Mathilda QuasiNewton BEATS scipy BFGS at n=1000 (1.65s vs 7.38s, 0.22x).
  * Mathilda LBFGSB SLOWER than scipy L-BFGS-B at large n (402 vs 12ms, 33x):
    the O(mn) solve is cheap, so per-solve COMPILE setup (objective+gradient to
    bytecode) dominates — same as bench63 A6/A7. Gradient-> barely helps
    (0.398->0.367s), confirming compile, not symbolic-diff, is the cost.

## Docs & verify — DONE
- [x] src/info.c docstring (LBFGSB method + aliases)
- [x] docs/spec/builtins/numerical-calculus.md + docs/spec/changelog/2026-08-10.md
- [x] make check-c99 PASS; valgrind clean (all lost = macOS runtime baseline noise)

## STATUS: COMPLETE. All suites green, no regression, C99-clean, valgrind-clean.
Not committed (awaiting user). M2 (exact GCP) deferred as documented follow-up.
