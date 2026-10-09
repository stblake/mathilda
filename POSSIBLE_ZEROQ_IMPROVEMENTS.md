# PossibleZeroQ / `zero_test` — improvement backlog

Follow-up items for the hybrid symbolic-numeric zero recogniser
(`src/zero_test.c`, `src/zero_test.h`; design in `ZERO_RECOGNISE_PLAN.md`).
Each entry: a minimal repro, the diagnosis, the observed impact, and a
concrete direction. These are **deferred** — logged here so a later session can
pick them up without re-deriving the analysis.

---

## Resolved & removed

Their full write-ups have been retired from this file; the fixes and their
regression tests (`tests/test_zero_test.c`) are the permanent record.

- **#1 — Gaussian × Erf residual spin.** RESOLVED 2026-09-07 (exp-combining
  Stage-0.5 normalisation; `exp_exponent_is_nonlinear`,
  `zt_normalize_exp_kernels`; test Group 16; commit `45bab708`).
- **#3 — `PossibleZeroQ[-E^(-a x)]` false positive.** RESOLVED 2026-09-26,
  v0.204 (structural non-zero certificate `provably_nonzero`, Stage 0b; test
  Group 17; commit `ce5463dc`).
- **#4 — special-function magnitude-decline residue (`PolyLog[2, x^2]+1`).**
  RESOLVED 2026-10-09, v0.323 (Stage-3 re-draw broadened to
  `residue_is_numeric_decline`; test Group 19; commit `b2f34965`).
- **#2 (PossibleZeroQ subset) — overflow-abort wrong `True`.** RESOLVED
  2026-09-26, v0.208 (IEEE-overflow re-draw in `evaluate_rung` /
  `sz_trial_shelled`; test Group 18; commit `7d87e186`). A residual *`Simplify`*
  deficiency from the same family remains open — see #2 below.
- **#6 — two-symbol coupling inequality unsound `False`.** RESOLVED
  2026-10-09, v0.325. `PossibleZeroQ[Max[x, y] - x, Assumptions -> x > y]`
  returned `False` though it is identically `0` on `{x > y}`. The per-symbol
  sampler drops a fact coupling two symbols, and `fact_keeps_false_sound`
  wrongly kept `False` for every inequality (the "full-measure region" rationale
  fails for a *piecewise* head like `Max`, a different function off the region).
  Fix: region-conforming **rejection sampling** — a coupling dense relation
  (`>`,`>=`,`<`,`<=`,`Inequality`,`Unequal`) over ≥2 of the expression's free
  symbols is collected into `conform_facts`, and `sz_trial` re-draws any
  assignment outside the feasible region (`assignment_conforms`, bounded by
  `ZT_CONFORM_MAX_TRIES`) before testing the point. A coupling *equality*
  (measure zero) keeps the `False → Unknown` downgrade. Engaged only when such a
  coupling is present, so all other draw streams are byte-for-byte unchanged.
  Test group `test_pzq_coupling_relations` in
  `tests/test_possiblezeroq_assumptions.c`.
- **#5 — huge-scale deep-cancellation screen false-negative.** RESOLVED
  2026-10-09, v0.324. At operand scale >= 2^`ZT_HUGE_SCALE_BITS` a
  machine-precision "obvious non-zero" is no longer trusted directly (a
  `(1 - near-1)` subtraction inside Tanh/Cosh/… loses >52 bits, so a genuine
  cancellation ZERO leaves a large rung-0 residual indistinguishable at machine
  precision from a true non-zero): the point is sent to the full MPFR ladder
  (`decide_numeric_huge_scale`), which trusts the residual's shrink/plateau
  trend — the cancellation collapses to True, while a genuine non-zero plateaus
  and stays False. `screen_point` defers huge-scale points to that ladder rather
  than settling a spurious whole-test False. Test Group 20.

---

## 2. `Simplify` hangs on a sum of exponentials with widely-separated real rates

**Status:** the `PossibleZeroQ` subset is RESOLVED (v0.208, see above); a
residual **`Simplify`-search** deficiency remains **OPEN**. Note this open item
is a `Simplify` bug, *not* a `PossibleZeroQ` one — `PossibleZeroQ` must never
call `Simplify`.

**Minimal repro:** `Simplify[Exp[-10 t] + Exp[-100 t]]` does not return (an 8 s
`TimeConstrained` aborts it); `Simplify[Exp[-t] + Exp[-2 t]]` returns instantly.

**Diagnosis.** `Simplify` invokes an equivalence zero-test (`simp_search`) that
numericalises the two terms at sample points. `E^(-10 t)` against `E^(-100 t)`
spans a huge dynamic range at a generic `t` (at `t=1.1`, `e^{-11}` vs `e^{-110}`
— a ratio of ~10^43), so the probe sees `small ± tiny` and the precision ladder
climbs to its ceiling without deciding. The separation, not the magnitude, is
the trigger: rates within ~1 order (-1,-2) never provoke it.

**Impact.** Constant-coefficient linear ODE *systems* with a real but spread
spectrum — e.g. `x'=-50x+20y, y'=100x-60y` (eigenvalues -10, -100) — produce a
fundamental-matrix body that is exactly such a sum. **Worked around (M27):**
`dsolve_linsys_tidy` (`src/calculus/dsolve_linsys.c`) routes any
exponential-carrying body through `Expand` rather than `Simplify`; the result is
back-substitution-verified regardless, so the cosmetic loss is harmless. A
direct `Simplify` of a wide-spectrum exponential sum still hangs.

**Suggested direction (core fix, deferred).** Before the numeric ladder, factor
out the dominant exponential (`E^(-10 t) + E^(-100 t) = E^(-100 t)(E^(90 t) + 1)`)
or compare terms in log-magnitude so a term negligible at the sample is not
differenced against a large one; or cap the ladder / return `UNKNOWN` when
per-point dynamic range exceeds the working precision.

### Cross-references (#2)

- `src/calculus/dsolve_linsys.c` — `dsolve_linsys_tidy` (the `Expand`-not-
  `Simplify` workaround and its rationale comment).
- DSolve §2.2.7 corpus systems 636, 650 (`DSolve_test_status/DE_examples_227.m`).
