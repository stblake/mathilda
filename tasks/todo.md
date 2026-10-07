# Mellin-transform Integrate: close 7 cases + fix ArcTan[a/x] sign bug (v0.300)

All edits in `src/calculus/integrate_ramanujan.c`; tests in `tests/test_integrate_ramanujan.c`.

## Implementation
- [ ] 1. Sign-bug fix: monomial Jacobian `1/k` → `1/|k|` (line ~1116)
- [ ] 2. `rec_expintegrale` — ExpIntegralE[n,a x] → a^-s Γ(s)/(s+n-1), strip Re s>0 ∧ Re(s+n)>1
- [ ] 3. `rec_trigpow` — Sin/Cos[a x]^k linearization (reuse sinpow_term w/ symbolic sv)
- [ ] 4. `conv_exp_trig` — Exp[-a x]{Sin,Cos}[b x] in rec_convolution
- [ ] 5. Coth[a x]-1 reduce rule (+ "Coth" guard) reusing rec_expgeom
- [ ] 6. `rec_arctan_sq` — ArcTan[a x]^2 → PolyGamma closed form
- [ ] 7. Register new recognizers in try_recognizers / rec_convolution

## Tests
- [ ] ArcTan[a/x] corrected (positive); regression: Sin[ax^2],Cos[ax^2],ArcCot,Log[1+a^2x^2]
- [ ] New: Sin^2, Sin^3, Exp·Sin, Exp·Cos, Coth-1, ArcTan^2, ExpIntegralE (2 cases)

## Verify / conventions
- [x] build + re-run all 12 pasted inputs — all correct; ArcTan[a/x] now POSITIVE
- [x] integrate_ramanujan_tests green (+ test_mellin_stress_set); definite suites green
      (intrep/diffunderint FAILs are PRE-EXISTING on clean tree, confirmed via stash)
- [x] independent oracle: 9/12 NIntegrate diff=0; 5,6,8 exact vs Pi/2, 0.794569..., 2/13
- [x] valgrind: identical leak profile to pre-existing Mellin baseline; no new-fn frames
- [x] version.h 0.299→0.300  (tag v0.300 pending commit)
- [x] docs/spec/builtins/calculus.md (table + prose) + changelog 2026-10-05.md
- [x] check-c99, check-messages clean; graph rebuilt

## Review
All 7 changes landed in src/calculus/integrate_ramanujan.c (general mechanisms,
no per-example patches). The ArcTan[a/x] "reference" in the paste was itself the
buggy (negated) value; the test asserts the corrected positive form. ExpIntegralE
returns a ConditionalExpression under the pasted (insufficient) assumptions — the
honest result (Re(s+n)>1 not provable), matching Mathematica. NOT yet committed
(awaiting user go-ahead for commit + git tag v0.300).
