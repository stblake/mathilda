# PossibleZeroQ — Round of Improvements (2026-07-02)

## Key finding (confirmed with the user)
B1 (`D[Integrate[f,t,"GoursatAlgebraic"],t] - f`) returning `False` is **CORRECT**, not a
bug: the antiderivative is only valid on the positive real axis, and the residual is
genuinely non-zero off that branch (`t=-3.4 → 0.20`). B1 shares the exact numeric
signature of `Sqrt[x^2]-x` and `Log[x^2]-2Log[x]`, both of which the suite *requires* to
be `False`. No sampling policy can flip B1 to `True` without breaking those. → NOT a PZQ
defect; keep `diff_back_ok`'s positive-real domain-aware gate in integrate_goursat.c.

## Scope (agreed)
- [x] Confirm B1 is not a bug (numeric evidence)
- [ ] **B2 performance** — early-stop the precision ladder once the residual has shrunk
      far below any deep-cancellation floor. Confirm phase = the whole cost (~0.8s/sample
      climbing to 1000 bits on a 60k-leaf tree). Target: 4.3s → <1s, no correctness loss.
- [ ] **Zippel-based rigor** — for *pure rational functions* (only free symbols, rationals,
      Plus/Times/Power[.,int]), a normalized non-zero numerator is a rigorous `False`.
- [ ] **Robustness + extensive stress tests**.
- [ ] **Docs/framing** — reclassify B1; note B2 fixed; changelog + docs/spec update.

## Out of scope (agreed)
- Genuine algebraic-function-field zero testing — deferred.

## Review (done 2026-07-02)
- [x] **B1 confirmed correct** (not a bug); docs reclassified, `diff_back_ok` retained.
- [x] **B2 perf**: deep-zero ladder early-exit (`ZT_DEEP_ZERO_BITS`) + reuse the
      precision-independent operand scale across rungs. 4.8s → 1.99s, verdicts unchanged.
- [x] **Zippel rigor**: `is_pure_rational_function` gate → decide_rational returns a
      rigorous `False` for pure rational functions (exact, deterministic, no sampling).
      `zero_test_decide` now trusts that `False`. Named algebraic constants excluded.
- [x] **Tests**: +21 cases (groups 12–14) in test_zero_test.c; full suite green in 1.3s.
- [x] **No regressions**: goursat/jeffrey/intrat/radical_simplify/simplify(143)/
      integrate_dispatch/chebychev all pass. solve_tests denest fail is PRE-EXISTING
      (identical on baseline build).
- [x] **valgrind**: PZQ paths identical to `Sin[1.0]` baseline (13,440B/420 blocks, all
      dyld/Accelerate init noise — zero Mathilda-src leaks).
- Docs: expression-information.md (stages 2–3), changelog 2026-06-29.md,
      POSSIBLE_ZEROQ_FAILURES.md.

## Changes (src/zero_test.c)
- `#define ZT_DEEP_ZERO_BITS 96` + early-return in the ladder shrink branch.
- `evaluate_rung` gains an `in_scale` param; higher rungs reuse the machine-precision scale.
- `is_pure_rational_function()`; `decide_rational` → rigorous FALSE; `zero_test_decide`
  trusts non-UNKNOWN from Stage 1.
