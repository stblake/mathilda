# DSolve M26 — §2.2.6 corpus (Problems 501–600) + general-forcing / DiracDelta

Plan file: `/Users/user/.claude/plans/let-s-continue-our-implementation-stateless-biscuit.md`

## Baseline (measured 2026-09-09)
`DE_examples_226.m`: 100 records (74 scalar / 47 IVP + 26 systems).
Corpus run: **57 PASS, 26 SKIP, 17 UNEVAL, 0 FAIL**.
The 17 U: forcing family 561–575 (15), the 555 hang, and 524 (Bessel verify quirk).

## Scoring insight
Branch verdict UNK → PASS (only BAD→FAIL, UNFIT→UNEVAL block). So arbitrary-`f`
scores PASS iff the **IC fit succeeds** (no leaked C[k]); DiracDelta scores PASS on
trust, so the in-method probe verify is the correctness guarantee.

## Corpus / converter  [done]
- [x] curl §2.2.6 HTML (WebFetch 403s); parse — 100 problems.
- [x] Converter `\delta(arg)→DiracDelta[arg]` fix (paren-guarded; bare `delta`
      Heun parameter preserved — verified via `convert_side`).
- [x] Generate `DSolve_test_status/DE_examples_226.m`.
- [x] Baseline measured.

## Solver work (ordered; each lands 0 FAIL)
- [x] 1. `integrate.c` — equal-limits rule `Integrate[_,{s,a,a}]→0`.
- [x] 2. `integrate_dirac.c`/`.h` (new) — DiracDelta sifting; called first in real-axis
      branch. (makefile auto-discovers src/calculus/*.c; added to tests COMMON_SRC.)
- [x] 3. `deriv.c` — variable-limit Leibniz rule + `HeavisideTheta' = DiracDelta`.
- [x] 4. `dsolve_common.c` — `dsolve_variation_of_parameters` definite-convolution
      fallback (fresh dummy `DSolve`vpS`; TrigReduce+Expand; undefined-fn gate).
- [x] 5. `distributions.m` (new, loaded from init.m) — H/DiracDelta value rules;
      `dsolve_constcoeff.c` DiracDelta-forcing gate. (Chose value rules over an
      in-method probe: the corpus verifier + probe-based unit tests are the gate.)
- [x] 6. `dsolve_common.c` — `dsolve_verify_body` spin guard (keep distributional residual).
- [x] 6b. `integrate.c` — skip improper/parametric methods on undefined-fn integrand
      (the real anti-hang fix: exp convolution 6s→0.3s per eval).
- [~] 7. 555 guard — reverted. The hang is upstream in `dsolve_linear1`'s solve of the
      regular-singular reduced eqn (pre-existing; U in baseline). Left as a documented
      residue (corpus 20s timeout handles it). ExactODE/spin-guard Ei attempts didn't fire.

## Wiring / docs
- [x] `tests/CMakeLists.txt` — `dsolve_corpus_2_2_6_tests` (gate baseline 2).
- [x] `tests/test_dsolve.c` — `t_m26_distributions/_impulse_forcing/_general_forcing`.
- [x] `DSolve_test_status/STATUS.md` §2.2.6 block; `README.md` row; `reports/2.2.6.{md,tsv}`.
- [x] `DSOLVE_PLAN.md` M26 entry.
- [x] `docs/spec/builtins/calculus.md` (D Leibniz, Integrate DiracDelta sift, DSolve forcing)
      + changelog `docs/spec/changelog/2026-09-07.md`.

## Verification
- [x] §2.2.6: **57→72/74 scalar, 0 FAIL** (2 residues: 524 slow-Bessel, 555 singular).
- [x] 561 convolution, 564 `½ Sin[2t] HeavisideTheta[t]`, 569 resonant, systems — all correct.
- [x] No regression: dsolve_corpus 2.2.1–2.2.5 gates held; series/reduce/integrate
      (dispatch/diffunderint/ramanujan/symmetry/newton_leibniz)/trigreduce; dsolve_tests +
      dsolve_stress_tests all pass.
- [x] `make check-c99` clean; `dsolve_corpus_2_2_6_tests` ctest **Passed** (gate baseline 2).
- [x] valgrind: new code (vp_definite_convolution / integrate_dirac / Leibniz) audited
      leak-free and appears in NO leak stack. A small per-call leak exists in the
      forced-equation cascade (`dsolve_factorable_try` / `poly_content` via ds_subst) —
      **pre-existing** (a pre-M26 forced solve like 534 leaks the same way), out of M26 scope.
      (macOS valgrind baseline ~13.4 KB noise; Linux CI is the definitive check.)

## Review
**M26 delivered: §2.2.6 corpus (Problems 501–600) added; 57→72/74 scalar solved, 0 FAIL.**
The whole general-forcing / DiracDelta family (561–575) now solves via Green's-function
variation of parameters (definite Duhamel convolution + DiracDelta sifting under Integrate
+ H/DiracDelta value rules + D Leibniz). Two pre-existing residues (524 slow Bessel, 555
singular-reduction series) remain, both UNEVAL with no wrong answers. Not committed
(awaiting user).
