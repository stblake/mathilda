# Task: Seamless Interval[...] support for special & piecewise functions

## Plan
Make special functions (Erfi, ExpIntegralEi, PolyLog, ...) and piecewise
functions (UnitStep, Ramp, Round, IntegerPart) thread over Interval[...] the way
elementary functions already do — rigorous enclosures, symbolic fallback.

## Steps
- [ ] `src/interval.c`: general derivative-sign certifier `interval_thread_call`
      (+ helpers iv_subst / iv_eval_at / iv_range_sign, depth guard).
- [ ] `src/interval.c`: bespoke rows in `interval_apply_function` for Erfi,
      InverseErf, InverseErfc, ProductLog, HarmonicNumber, UnitStep, Ramp,
      Round, IntegerPart.
- [ ] `src/interval.h`: declare `interval_thread_call`.
- [ ] `src/eval.c`: central hook after builtin returns NULL (NumericFunction +
      interval arg -> interval_thread_call).
- [ ] `src/piecewise.c`: update the "deliberately not threaded" comment.
- [ ] `tests/test_interval.c`: new coverage cases.
- [ ] Docs: docs/spec/builtins interval reference + changelog 2026-08-24.
- [ ] Verify: build, re-probe, `make check-interval`, full test suite,
      confirm oscillatory heads stay symbolic and terminate.

## Review

Implemented. Special & piecewise functions now thread over `Interval[...]`.

**Mechanism (elegant + rigorous):**
- `interval_thread_call` (src/interval.c) — general derivative-sign certifier:
  interval-evaluates `D[f[x],x]` over each pair; entirely ≥0 → increasing,
  ≤0 → decreasing; else symbolic. Reuses `D[]` + interval evaluator, no
  per-function monotonicity analysis. Wired into src/eval.c step "5b" (after
  builtin NULL, gated on ATTR_NUMERICFUNCTION + interval-arg scan).
- Depth cap `IV_CERTIFY_MAX_DEPTH`=4 + counter-raised-across-whole-section:
  chains (PolyLog[3]) resolve; oscillatory chains (Bessel) stay symbolic + fast.
- `iv_range_sign` only sign-tests a threaded Interval or concrete number (never
  numericalizes a symbolic derivative → no recursion, fast fallback).
- Bespoke rows: InverseErf, InverseErfc, ProductLog, HarmonicNumber, Erfi;
  piecewise UnitStep/Ramp/Round/IntegerPart (non-decreasing).

**Rigor verified:**
- `make check-interval` — containment 23300→29970 (new heads added), all pass.
- Interior-sampling spot-checks: every enclosure contains all sampled values.
- Discontinuous non-monotone (Mod, FractionalPart) correctly stay symbolic
  (D leaves them as unevaluated Derivative forms).
- valgrind: definitely-lost identical at 1 vs 40 iters (fixed macOS objc/dyld
  startup noise), zero interval.c/eval.c frames → certifier path leak-free.
- 13 targeted test binaries pass (interval, eval, core_algebra, + each touched
  head); full 454-binary suite run for regressions.

**Deliberately symbolic (documented):** oscillatory heads (Bessel*, Sinc,
SinIntegral, Fresnel*, Airy*) and sawtooth FractionalPart — no certifiable
monotone bound; safe (never wrong).
