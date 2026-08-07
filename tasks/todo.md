# Assumption-aware Limit — todo

Plan: /Users/user/.claude/plans/in-the-same-way-encapsulated-minsky.md

## Phase 0 — Stress corpus (fitness function) ✅
- [x] New tests/test_limit_stress.c (62 cases, 11 categories)
- [x] Extend tests/test_limit_assumptions.c (2 wiring functions)
- [x] Wire tests/CMakeLists.txt (limit_stress_tests)
- [x] Build + confirm [C] pass, [N] FAIL (46 failing: 41 stress + 5 assumptions)

## Phase 1 — Wiring + generalized early dispatch ✅
- [x] #include simp.h; add const AssumeCtx* assume to LimitCtx (+ 13 brace ctors)
- [x] limit_effective_assumptions (option + ambient; Automatic-as-none; inconsistent->NULL)
- [x] single-exit assume_ctx_free (goto cleanup); thread into run_iterated
- [x] de-static read_dollar_assumptions (simp_builtins.c + simp.h)
- [x] limit_real_base_power via assume_known_gt/lt(ctx, base, 1)

## Phase 2 — Sign oracle ✅
- [x] assume_known_gt/lt/gt_expr in simp_assume.c + simp.h (+ Inequality decomposition)
- [x] literal_sign_ctx in limit.c; ctx Log[b] arm; exponent_sign_ctx difference arm
- [x] rewire read_leading_term_limit + layer3_rational

## Phase 3 — Growth exponent + monomial handler ✅
- [x] layer_param_power (x^exp monomial) as a compute_limit LAYER (runs first)
- [x] layer_param_monomial_substitute (t = x^a) for rational functions of x^a
- [note] growth_exponent_upper threading NOT needed — bounded envelope resolves
  x^(-p) via layer_param_power on the recursive compute_limit

## Phase 4 — Compose-at-infinity direction ✅
- [x] falls out of Phase 2: sub-limit c x returns signed Infinity, ArcTan/Tanh fold

## Phase 5 — Domain -> sign polish ✅
- [x] PositiveIntegers/NonnegativeIntegers/NegativeIntegers/NonpositiveIntegers => sign

## Docs & verification ✅
- [x] docs/spec/builtins/calculus.md Limit Assumptions bullet + examples
- [x] docs/spec/changelog/2026-08-03.md feature summary
- [x] Limit docstring refresh (in limit.c, verified via ?Limit)
- [x] make check-c99 clean; all suites 0 FAIL; valgrind: no leak stack touches new code

## Review
- Corpus: 62-case tests/test_limit_stress.c (11 categories) + 2 wiring functions in
  test_limit_assumptions.c. Started 46 FAIL -> 0 FAIL.
- Development landed as: shared simp_assume additions (threshold predicates,
  chained-inequality decomposition, domain->sign); a borrowed AssumeCtx threaded
  through LimitCtx; literal_sign_ctx / exponent_sign_ctx; two new cascade layers
  (layer_param_power, layer_param_monomial_substitute).
- Regressions: limit/gruntz/oscillatory/possiblezeroq(x2)/simplify/element/series/
  powerexpand/fullsimplify(x2)/nlimit/integrate(newton-leibniz,dispatch) all 0 FAIL.
- Sound-only: NULL ctx == legacy path byte-for-byte; verdicts only on entailment.
- Not commited (awaiting user).

## Series honours Assumptions / Assuming[] (2026-08-07)
- Goal: replicate the Limit-Assumptions work for Series + a 100+ case stress corpus.
- Done: `series_effective_assumptions` (HoldAll → eval option first) + borrowed
  `const AssumeCtx*` threaded through SeriesCtx/do_series_single/series_expand;
  sign-of-x from ctx (dropped assumption_sign_of/lit_real_sign + int x_sign);
  non-analytic kernels (Abs/Sign/UnitStep/Conjugate of x, and Sqrt[g^2]→±g);
  final-coefficient cleanup via apply_assumption_rules (published in simp.h; added
  (x^2)^r→x^(2r) both signs, Sign/Conjugate rules).
- Corpus: tests/test_series_stress.c (110 cases, 10 cats) + test_series_assumptions.c
  (wiring) — every expected verified against the binary. All 0 FAIL.
- Regressions: series/series_twoterm/nseries/limit_stress/limit_assumptions/
  possiblezeroq_stress/assuming all 0 FAIL; simplify has 1 pre-existing cosmetic
  FAIL (x^2^(3/2) stale baseline, unrelated). check-c99 clean. NULL ctx == legacy
  byte-for-byte. No new leaks (valgrind: my funcs never allocate a lost block).
- Not committed (awaiting user).
