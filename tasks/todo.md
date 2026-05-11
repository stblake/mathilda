---
title: IntegrateRational closed-form gaps
date_started: 2026-05-11
status: in progress
---

# IntegrateRational closed-form gaps (2026-05-11)

Five failing `Integrate`IntegrateRational[...]` cases.  All five
reproduce in picocas REPL; all five integrate cleanly via the
Mathematica .m source (verified with wolframscript), so the gaps are
in the C port at `src/intrat.c`, not in the underlying algorithm.

## Root-cause map

| # | Input | Root cause |
|---|-------|------------|
| 1 | `1/(b - a x^3)` | Cubic Q-in-t.  `logtoreal_dispatch` has no deg-3 branch; falls to RootSum. |
| 2 | `1/(x^5+1)` | Palindromic-quartic Q-in-t `5(t^4-t^3+t^2-t+1)`.  `logtoreal_dispatch` deg-4 only handles biquadratic & Sophie-Germain.  Falls to RootSum. |
| 3 | `x^4/(a+b x^3)^2` | Same as #1 — cubic Q after Hermite reduction. |
| 4 | `(-1+x^2)/(1-2x^2+2x^4)` | Biquadratic Q-in-t `2t^4-2t^2+1` with negative inner discriminant.  Current deg-4 c1z&&c3z branch bails on neg disc; falls to `expand_simple_rootsum`, which produces complex form. |
| 5 | `1/x/(1-2x^2+2x^4)` | Post-processing: `c * piece_int` stays held as `Times[c, Plus[…]]`; also `Log[c p]` not stripped when `FreeQ[c, x]`. |

## Phases (independently testable)

### A — Output post-processing (#5)
- Distribute scalar Times-over-Plus in the final accumulator.
- Strip `Log[c · p] -> Log[p]` when `FreeQ[c, x]`.

### B — Negative-inner-disc biquadratic (#4)
- Extend `logtoreal_dispatch` deg-4 c1z&&c3z branch: when inner disc
  `c2^2 − 4 c0 c4 < 0` (provably under positive-symbol assumption)
  and `c0 c4 > 0`, factor over R as
    `c4 (t^2 + α t + Sqrt[c0/c4])(t^2 − α t + Sqrt[c0/c4])`
  with `α = Sqrt[2 Sqrt[c0/c4] − c2/c4]`.
- Dispatch each factor through `logtoreal_quadratic`.

### C — Cubic (#1, #3)
- Add deg-3 branch to `logtoreal_dispatch`.
- C.1: nth-root cubic (c1 = c2 = 0).  Real root
  `r0 = −(c0/c3)^(1/3)`; factor `c3 (t − r0)(t^2 + r0 t + r0^2)`;
  dispatch the linear and quadratic factors.

### D — Palindromic quartic (#2)
- Detect `c0 t^4 + c1 t^3 + c2 t^2 + c1 t + c0`.
- `u = t + 1/t` ⇒ `c0 u^2 + c1 u + (c2 − 2 c0) = 0`.  Factor as
  `c0 (t^2 − u_+ t + 1)(t^2 − u_− t + 1)` and dispatch through
  `logtoreal_quadratic`.

### E — Tests
- Add picocas tests for all 5 cases (no `RootSum`/`Function` head in
  output; `D[result, x] − integrand` simplifies to 0).
- Extend `IntegrateRationalTests.m` for Mathematica parity.

## Review

All 5 phases shipped.

### Phase A — Output post-processing (#5) ✓
Added `intrat_distribute_plus(e)` (uses `expr_expand`) and
`intrat_strip_log_constants(e, x)`.  Both run before AND after
`intrat_log_to_arctanh` in `intrat_integrate_rational`.

### Phase B — Negative-inner-disc biquadratic (#4) ✓
Extended `logtoreal_dispatch` deg-4 `c1z && c3z` branch with a
Sophie-Germain-with-c2 sub-branch that factors as
`c4 (t^2 + α t + β)(t^2 − α t + β)` with `β = Sqrt[c0/c4]`,
`α = Sqrt[2 β − c2/c4]`.

### Phase C — Cubic (#1, #3) ✓
New deg-3 branch in `logtoreal_dispatch` for nth-root cubics
`c3 t^3 + c0`.  Real-root selection under positive-symbol or
numeric-sign assumption; dispatch linear + quadratic factors.

### Phase D — Palindromic quartic in logtoreal_dispatch ✓ (restricted)
Detect scaled-palindromic `c4 c1^2 == c3^2 c0`, but only fire when
r = 1 (pure palindromic).  Scaled cases would need a smarter
`LogToAtan` (current C port hangs on nested-radical S).

### Phase D2 — Palindromic quartic in NaiveLogPart (#2) ✓
New `expand_palindromic_quartic_real(a_t, d_t, dd_t, t, x)` builds
the real form directly from `(a(α) / d'(α)) · Log[x − α] + …` per
conjugate root pair.  Bypasses LogToReal / LogToAtan entirely; the
nested-radical S that overwhelms LogToAtan on the LRT path never
gets constructed here.

### Phase E — Tests ✓
`tests/test_intrat.c` adds five new regression tests (one per
issue) using `assert_closed_real` (no `RootSum` / `Function`
leak) + `assert_integral_numeric_ok` for nested-radical cases.
All five pass.  `IntegrateRationalTests.m` gains five new
{integrand, var, optimal-antiderivative} entries.

### New helpers added
- `intrat_distribute_plus(e)` — Plus-over-Times distribution.
- `intrat_strip_log_constants(e, x)` — strip x-free constants out of Log args.
- `intrat_numeric_sign(e)` — N[]-based sign decision, fallback for
  symbolic positive-symbol bails on `Sqrt[5] − 5` style inputs.
- `expand_palindromic_quartic_real(a_t, d_t, dd_t, t, x)` —
  real-form expander for palindromic-quartic d in NaiveLogPart.

### Verification
- Build clean under `-std=c99 -Wall -Wextra -O3`.
- All five regressions return `Plus`-headed real-elementary form.
- Differentiation check passes (symbolic for #5, numeric for #1–#4
  whose nested radicals don't reduce by Cancel/Together).
- Full `tests/` suite: only pre-existing `logexp` + `poly` ordering
  failures (introduced by 2026-05-11's canonical comparator change),
  no new regressions from the closed-form work.

### Known follow-ups
- Cardano's formula for non-depressed cubics in `logtoreal_dispatch`.
- Scaled-palindromic dispatch (r != 1) — requires a cheaper LogToAtan
  for nested-radical coefficients, or a different routing strategy.
- The closed-form outputs aren't always Mathematica-minimal; the
  underlying form is correct but Sqrt simplifications like
  `Sqrt[1/16 (1 + Sqrt[2])] -> Sqrt[1+Sqrt[2]]/4` are left to a
  smarter `intrat_simp_pos_sqrt` (or full `Simplify`).

## 2026-05-11 rev 2 — generalised nth-root coverage

User flagged that the cubic-specific fix above was too narrow.
Adjacent inputs all leaked `RootSum`:

  - `1/(b - a x^4)`, `1/(b - a x^5)`,
    `1/(b - a x^6)`, `1/(b + a x^6)`,
    `1/(b + a x^8)`, `1/(b + a x^9)`.

### Generalisation
Replaced the per-degree branches in `logtoreal_dispatch` with one
helper `logtoreal_nthroot_sparse(base, deg, s, x, t)` that closes
`c_n t^n + c_0` (all intermediate coefficients zero) for every
n ≥ 3 via the standard cyclotomic factorisation over R.  Subsumes
the previous deg-3 nth-root and deg-4 Sophie-Germain branches.

Algorithm: given `q = -c_0/c_n` and `r = |q|^(1/n)`,
- q > 0: real roots ±r (r always; -r iff n even), conjugate pairs
  for k = 1..⌊(n-1)/2⌋ at angles 2πk/n.
- q < 0: real root -r iff n odd, conjugate pairs at angles
  (2k+1)π/n.
Each conjugate pair contributes a real quadratic
`t^2 - 2 r cos(θ) t + r^2` with discriminant -4 r^2 sin²θ < 0, so
`logtoreal_quadratic`'s ArcTan branch always fires.

### Tests
Six new regressions in `tests/test_intrat.c`:
`test_closed_nth_root_quartic_minus`, `_quintic_minus`,
`_sextic_minus`, `_sextic_plus`, `_octic_plus`, `_nonic_plus`.
Each asserts no RootSum/Function head + numerical derivative
match at two sample points.  All pass.

### Regression status
Full test suite (105 binaries): 90 pass / 15 fail both before and
after this change.  Identical failure set (pre-existing canonical
ordering issues from commit cc8b164).  Zero new regressions.
