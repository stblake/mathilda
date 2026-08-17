# Solve[] review — stress corpus + test-driven fixes

Plan: `/Users/user/.claude/plans/let-s-do-a-review-indexed-pond.md`

## Phase 1 — Stress-test corpus + back-substitution harness  [DONE]
- [x] `tests/solve_check_prelude.m` — `solveVerdict`/`solveCheckCode`/`solveReport` (HoldAll; count + back-sub + domain; Root/ConditionalExpression/Modulus shape-aware). Thread-based sampling avoids Function-closure gap.
- [x] `tests/solve_corpus.m` — 62 case records, categories A–H, refNotes
- [x] `tests/solve_corpus_run.m` — standalone driver
- [x] `tests/test_solve_corpus.c` — fork-per-case C runner cloned from test_intrat_corpus.c, baseline gate
- [x] `tests/CMakeLists.txt` — register solve_corpus_tests
- [x] `SOLVE_FAIL_BASELINE = 10` (D1 Modulus×4, D2 kernel×2, D3 Rationals×2, D4 multivar×2). 52/62 pass.

## Phase 2 — D3 Rationals domain (tiny)  [DONE]
- [x] solvepoly.c domain switch (rationals_only) + is_rational_like tail filter
- [x] updated test_solve.c: new test_rationals_domain; Algebraics stays unevaluated
- [x] corpus baseline 10 -> 8; solve_tests PASS. (solvenlsys Rationals deferred to Phase 6 alongside multivar)

## Phase 3 — D6 mixed-trig simplification (tiny)  [DONE]
- [x] solvetrig.c aggregate tail: Simplify RHS only when it removes a complex Log (descend ConditionalExpression), leaf-count guard. tan1 cleaned to Pi/4+Pi C[1]; clean sums not churned.
- [~] cos-sin (Cos==Sin) still shows complex-log: Simplify/Arg can't reduce Log of unit-modulus complex; that's a Simplify subsystem gap, out of scope.

## Phase 4 — D1 Modulus (critical bug)  [DONE]
- [x] src/solvemod.{c,h} (residue enumeration, refuse systems/non-poly/out-of-range); Solve`SolveModular builtin
- [x] SolveOpts.modulus; apply_option wires Modulus; pre-pass in solve.c (goto solve_finish tail); solvemod_init
- [x] tests: COMMON_SRC += solvemod.c; new test_modulus_domain (4 solve + 2 refuse); Modulus docstring
- [x] corpus baseline 8 -> 4; solve_tests PASS

## Phase 5 — D5 VerifySolutions (thin PossibleZeroQ wrapper)  [DONE]
- [x] SolveOpts.verify_on; apply_option wires VerifySolutions->True
- [x] post-dispatch filter in solve.c: solution_verifies (zero_test_decide==FALSE drops; Root/param/ConditionalExpression kept) + verify_solutions_filter; leak-safe early-return
- [x] docstring updated (no longer "reserved"); test_verify_solutions (keeps poly/radical/Root); solve_tests + corpus PASS

## Phase 6 — D4 single-equation many-variables  [DONE]
- [x] try_single_eq_multivar in solve.c between linsys-NULL and nlsys (single Equal; earliest polynomial-solvable var; explicit rules only)
- [x] updated test_positive_dimensional to a genuine 2-eqn system; new test_single_eq_multivar; solve_tests + solvenlsys_tests PASS
- [x] corpus baseline 4 -> 2

## Phase 7 — D2 poly in one transcendental kernel  [DONE]
- [x] solvetrig_solve_poly_in_kernel: exp path (subst_exp_real_walk, E^(cx)->u^c) + generic-head path (collect_kernels, structural u-subst); solve_u_and_unwind reuses solvepoly + solveinv
- [x] wired in solve.c between solvetrig and solverad; prototype in solvetrig.h
- [x] new test_poly_in_kernel (Log^2, Log^3, decline E^x+x); corpus baseline 2 -> 0 (FULLY GREEN)

## Docs
- [ ] docs/spec/builtins/solutions-of-equations.md + docs/spec/changelog/2026-08-17.md
- [ ] Refresh Solve docstring (solve.c)
- [ ] make check-c99

## Docs  [DONE]
- [x] docs/spec/builtins/solutions-of-equations.md (dispatch list, Modulus/Rationals/VerifySolutions, Domains note)
- [x] docs/spec/changelog/2026-08-17.md (Solve review section)
- [x] Solve docstring in solve.c (options list + narrative); Modulus/VerifySolutions docstrings; Solve`SolveModular docstring
- [x] make check-c99 clean

## Review  [COMPLETE]
- All 62 corpus cases pass; SOLVE_FAIL_BASELINE ratcheted 10 -> 0 across phases.
- Full ctest suite: 219/219 pass, 0 regressions (incl. integration/sum/cherry that call Solve internally).
- valgrind: new code introduces ZERO leaks (baseline dyld/runtime noise identical with/without the new features).
- Files added: src/solvemod.{c,h}, tests/{solve_check_prelude.m, solve_corpus.m, solve_corpus_run.m, test_solve_corpus.c}.
- Files modified: src/solve.c (Modulus pre-pass, VerifySolutions filter, single-eq-multivar), src/poly/solvepoly.c (Rationals), src/solvetrig.{c,h} (poly-in-kernel + mixed-trig Simplify), tests/{test_solve.c, test_solvenlsys.c, CMakeLists.txt}, docs.
- Boundary respected: no inequalities/quantifiers/case-splits; Reduce territory untouched.
