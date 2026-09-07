# DSolve M17 — 2nd-order-linear special-function wave

## Part A — Affine → Gauss ₂F₁   ✅ + F-homotopy (exponent shift, needed for Gegenbauer/Jacobi/assoc-Legendre)
- [x] Extract canonical Gauss row into static `specialform_gauss_basis(Pc,Qc,xvar)`
- [x] Replace canonical Gauss block with a call to the helper (byte-identical output)
- [x] Add `hgc_num_ok(...)` numeric self-verify (in-segment sampling; instantiates x1,h params)
- [x] Add affine row: locate 2 finite RSPs (squarefree PolynomialLCM deg 2), map x=x1+h·s
- [x] **F-homotopy** `specialform_fhomotopy`: pull local exponents Y=s^r0(1-s)^r1 F → canonical Gauss
- [x] Decline deg L ≠ 2 / h==0 cleanly with full frees
- [x] Use Cancel[Together] (not Simplify) for Pt/Qt and pole extraction (Simplify drops a pole)
- [x] Legendre row: emit only ordinary (μ==0); associated (μ≠0) → affine 2F1

## Part B — Normal-form pre-pass   ✅
- [x] Extract 4 P==0 rows into `specialform_reduced_basis(Qc,xvar)`
- [x] Replace inline blocks with helper call (byte-identical output)
- [x] Add P≠0 normal-form pre-pass row: dsolve_normal_form → -r → reduced_basis → mu·base, gate sf_num_ok
- [x] Update DSolve`SpecialFunctionForm docstring

## Tests & docs
- [x] Unit tests in test_dsolve.c (t_m17_affine_gauss/gegenbauer/associated_legendre/normalform_bessel/declines)
- [x] New tests/test_dsolve_m17_stress.c + CMake add_executable/add_test (Gegenbauer/Jacobi/assocLeg/shifted/Bessel grids)
- [~] Reconcile corpus gate argv[3] 612 → measured N  (WAITING on corpus run)
- [x] DSOLVE_PLAN.md M17 entry + changelog 2026-09-07.md + calculus.md
- [~] STATUS.md scoreboard  (WAITING on corpus number)

## Gates
- [x] make -j build + make check-c99 green
- [x] dsolve_tests (incl 5 M17 units), dsolve_m5/m12/m14_stress, dsolve_stress, dsolve_m17_stress all green
- [~] ctest dsolve_corpus_2_1_2_tests: measuring N (running)
- [ ] valgrind leak-clean on new stress binary + dsolve_tests

## Side-finding (user asked to fix): Factor drops squared term  ✅ FIXED
- Factor[(5-3s)^2 - 30 s] = -5(-5+12s) WRONG (dropped 9s²); now = 9s²-60s+25.
- Root cause: flint_univariate_factor computed degree on raw (unexpanded) arg;
  get_degree_poly counts Power[base,k] only when base==var. Fix: expr_expand first.
- Also fixes Simplify (calls Factor). Regression test: test_factor_baseline test_unexpanded_power_input.

## Corpus regressions from the affine row / normal-form pre-pass  ✅ FIXED
- 11 apparent PASS→UNEVAL: 7 second-order latency (affine row slow on decline path),
  2 third-order $IterationLimit loops (900/1199), 1 pre-existing flaky (997, loops on HEAD too), 603 ok.
- Fixes: quad-size gate (80 leaves) in specialform_quad_roots; radical-RSP gate (sf_has_radical);
  sf_ct (Cancel[Together]) instead of Simplify in fhomotopy/gauss_basis; normal-form pre-pass
  size gate (negr≤50) + depth gate (g_dsolve_depth≤1, avoids OperatorFactor-recursion loop).
- Verified: all 10 real regressions solve <8s; 20+ improvements preserved; spherical-Bessel pre-pass kept.
