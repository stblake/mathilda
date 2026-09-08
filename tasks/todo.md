# PossibleZeroQ hang fix — exponential-combining normalisation

Fix the `PossibleZeroQ` hang on Gaussian × Erf residuals (item #1 of
`POSSIBLE_ZEROQ_IMPROVEMENTS.md`). Plan:
`/Users/user/.claude/plans/let-s-improve-the-recently-snug-pebble.md`.

## Core change — src/zero_test.c  [DONE]
- [x] Add `exp_exponent_is_nonlinear` + `expr_has_symbolic_exp_kernel` gate (narrowed to NON-LINEAR exponents).
- [x] Add `zt_normalize_exp_kernels(const Expr*)` → `expr_expand_all` or NULL.
- [x] Extract `zt_decide_core` from `zero_test_decide` body; wrap with normalisation.
- [x] Extract `zt_decide_assuming_core` from `zero_test_decide_assuming` ctx-body; wrap.
- [x] Update file/header docstrings.

## DSolve workaround removal — src/calculus/dsolve_common.c  [DONE]
- [x] Remove Erf/Erfi ExpandAll pre-pass (408-419).
- [x] Remove FALSE-path Power[E,ztexp] re-check (421-436) — §2.2.5 corpus + y'+xy==Exp[3x] still solve.

## Tests  [DONE]
- [x] test_zero_test.c: Group 16 "exponential-combining" (repro + preservation + must-be-False + stable).
- [x] New tests/test_possiblezeroq_expcombine_stress.c (ctest-registered #220, hard exit, 31 cases).
- [x] Register stress file in tests/CMakeLists.txt (no COMMON_SRC change).

## Key finding during implementation
- Blanket gate regressed a constant-coeff-ODE UC case (residual has incidental affine E^x that
  cancels; ExpandAll mangled the trig part into a sampler-hostile form → False). FIX: narrow gate
  to NON-LINEAR exponents (the tiny*huge overflow needs super-linear growth). Affine E^x untouched.
- F4/F5-style flat-product & D[b,x]-b false-positives are PRE-EXISTING sampler limits (ExpandAll
  is a no-op there); documented, not asserted.

## Docs  [DONE]
- [x] docs/spec/changelog/2026-09-07.md — PossibleZeroQ section (newest-first).
- [x] docs/spec/builtins/expression-information.md — Stage 0 normalisation note.
- [x] POSSIBLE_ZEROQ_IMPROVEMENTS.md — item #1 marked RESOLVED.

## Review
Fix landed as a value-preserving Stage 0.5 exponential-combining normalisation in
`src/zero_test.c`, applied at the top of both public entries (above the Stage-3
routing gates the failing input hit). Gate narrowed to NON-LINEAR exponents after a
blanket gate regressed a constant-coeff-ODE verify. DSolve workaround removed. New
Group 16 in test_zero_test.c + ctest-registered test_possiblezeroq_expcombine_stress.c
(31 cases). Broad regression sweep clean: zero_test, trigexp, PZQ assumptions/stress,
dsolve, dsolve_stress, dsolve_m14_stress, simplify (1 FAIL pre-existing, no exp, soft
assert), comparisons, integrate_risch_transcendental, integrate_fresnel, knowles_erf,
erf, erfi — all pass / unchanged. valgrind identical to baseline; check-c99 clean.
NOT committed (awaiting user).

## Verification  [DONE except docs]
- [x] Build main + tests (GCC clean, no unused-function warnings).
- [x] Repro fixed (PossibleZeroQ[res] fast True; was >20s hang / 7703 iter-limit floods).
- [x] zero-test/trigexp/assumptions/stress + DSolve/DSolve-stress suites: all pass, no verdict change.
- [x] DSolve cases 428, 482, y'+xy==Exp[3x] still solve.
- [x] valgrind: identical to baseline (0 new leaks/errors); make check-c99 clean.

## Review
(to be filled in)
