# Campaign: Solve 12000.org §2.1.2 ODE corpus in DSolve

Plan: `~/.claude/plans/let-s-devise-algorithmic-methods-graceful-magpie.md`
Corpus: Table 2.3, 1204 ODEs (Maple+Mma solve all; SymPy 157). Scalar-first,
measured-biggest-gaps-first, full campaign to parity.

## Phase 0 — Infrastructure & baseline (M15)
- [x] `tools/latex_ode_to_mathilda.py` — LaTeX→Mathilda converter (reads section HTML,
      array-preamble strip, function/indvar detection, subscripts, operatorname/textit).
- [x] Validated: all 1391 equations round-trip through Mathilda's parser (0 errors).
- [x] `DE_examples_2.m` (repo root) — 1204 records (1000 scalar + 204 systems), Get→1204.
- [x] `tests/test_dsolve_corpus.c` + `tests/dsolve_corpus_prelude.m` — fork-per-case
      harness, numeric back-subst self-verify, systems skipped, per-case TSV out.
- [x] Register `dsolve_corpus_tests` in `tests/CMakeLists.txt`; builds + smoke-tested.
- [x] `tools/dsolve_corpus_report.py` — buckets TSV by Maple class, ranks gaps.
- [x] Consolidated into `DSolve_test_status/` (corpora + harness + prelude + README +
      STATUS + reports/); CMake references `../DSolve_test_status/`, per-section add_tests.
- [x] Baseline: 385/1000 (38.5%), 0 FAIL, 6 crashes → ranked gaps in reports/2.1.2.md.
- [x] `DSolve_test_status/STATUS.md` — cross-session dashboard (2.1.2 numbers in).
- [x] **Crash hardening** (surfaced by corpus): (A) `dsolve_linear_normalize`
      rational-in-x gate (`ds_is_rational_in`) — kills GCD(-inf,1) stack overflow;
      298/876/879 now SOLVE. (B) `expr_to_mpolyq` zero-inverse guard — kills FLINT
      abort (269). Post-fix: 388/1000, gap 615→**612**, crashes 6→2 (both
      non-reproducing/flaky). All DSolve suites + check-c99 green.
- [x] DSOLVE_PLAN.md M15 entry + docs/spec changelog.
- [ ] Section **2.1.3**: user curls section3.html → convert → `DE_examples_3.m` →
      baseline → add `dsolve_corpus_2_1_3_tests`. (awaiting fetch)
- [ ] Chase intermittent crashes 2.1.2-208 / 2.1.2-983 if they become reproducible.

## Phase 1 — M16 ✅ DONE: 2nd-order-linear change-of-variable + recognizer
- [x] `cv_num_ok` instantiates free params → symbolic-degree Legendre-Cot solves
      (`dsolve_changevar.c`); guard `t_m16_legendre_symbolic`.
- [x] Pöschl-Teller / trig-potential recognizer → verifiable 2F1, in-method numeric
      self-verify (`dsolve_specialform.c`); guard `t_m16_poschl_teller`.
- [x] 388→**396 solved**, gap 612→604, 0 FAIL; all DSolve suites + check-c99 green.
      Baseline 605 (601 deterministic + margin for flaky fork crashes 208/872/983).
- [x] DSOLVE_PLAN M16 + STATUS + changelog + docstring updated.
Finding: change-of-*independent*-var substitutions (x^m/e^x/ln x/1/x) ~0 yield on
homogeneous-rational (need recognizers at symbolic exponents); canonical Gegenbauer/Jacobi
already Kovacic-solved.

## Next: M17 (measurement-ordered, not started)
- [ ] Abel Invariant Rational (64, deferred M13) — biggest single-method 1st-order gap
- [ ] 3rd/high-order operator factoring (110)
- [ ] 2nd-linear residue: Gegenbauer/Jacobi via affine→Gauss 2F1 (math validated),
      parabolic-cylinder, power→Bessel
- [ ] Section 2.1.3: awaiting user curl of section3.html
- [ ] Investigate flaky fork-harness crashes (208/872/983) if they become reproducible
Candidate targets (measured, ranked): Abel AIR (64, single coherent method, deferred
M13) vs 2nd-order-linear special-function families (222, fragmented: Heun /
parabolic-cylinder / power→Bessel / Pöschl-Teller) vs 3rd/high-order operator
factoring (87). Each: 3-fn contract + cascade wiring + 10-20 anti-overfit family +
DSOLVE_PLAN milestone + docs + re-run corpus (ratchet 612 down) + suites + check-c99.

## Phase 1..N — Method waves (M16+), measurement-ordered
- [ ] Wave A — Abel Invariant Rational (revive M13), `dsolve_abel_air.c`
- [ ] Wave B — orthogonal-poly recognizers (symbolic degree), extend `dsolve_specialform.c`
- [ ] Wave C — 2nd-order nonlinear μ integrating-factor reduction, `dsolve_mu_reduce.c`
- [ ] Wave D — higher-order operator factoring (Beke/2nd-order right factors)
- [ ] Wave E — solvable-for-y/x general first-order, `dsolve_solvefor.c`
- [ ] Wave F — 2nd-order linear change-of-variable extensions (t=x^k, e^x, ln x, Möbius)
- [ ] Wave G+ — long tail (Emden–Fowler, elliptic, Liénard, rational/NONE)

## Notes
- Each wave: 3-fn contract + cascade wiring + anti-overfit stress family + DSOLVE_PLAN
  milestone + docs/changelog + re-run corpus (ratchet baseline) + ctest + check-c99 + valgrind.
- 204 systems recorded now, deferred to a later campaign.
