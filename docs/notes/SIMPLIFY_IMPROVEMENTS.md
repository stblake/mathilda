# Simplify — improvement backlog

Follow-up items for the general `Simplify` pipeline (`src/simp/simp_search.c`).
Each entry: a minimal repro, the diagnosis, the observed impact, and a concrete
direction. These are **deferred** — logged here so a later session can pick them
up without re-deriving the analysis.

> **Relocated 2026-10-16** from the retired `POSSIBLE_ZEROQ_IMPROVEMENTS.md`
> (item #2). The `PossibleZeroQ` subset of that item — the overflow-abort wrong
> `True` — was resolved separately (v0.208; regression `tests/test_zero_test.c`
> Group 18). What remains below is a **`Simplify`-search** deficiency, *not* a
> `PossibleZeroQ` one: `PossibleZeroQ` must never call `Simplify`.

---

## 1. `Simplify` hangs on a sum of exponentials with widely-separated real rates

**Status:** OPEN. Re-verified still reproducing 2026-10-16, v0.420.

**Minimal repro:** `Simplify[Exp[-10 t] + Exp[-100 t]]` does not return (an 8 s
`TimeConstrained` aborts it → `$Aborted`); `Simplify[Exp[-t] + Exp[-2 t]]`
returns instantly as `(1 + E^t) E^(-2 t)`.

**Diagnosis.** Two accounts, because profiling the current build complicates the
original one:

- *Original (dynamic-range) account.* `Simplify` invokes an equivalence
  zero-test that numericalises the two terms at sample points. `E^(-10 t)`
  against `E^(-100 t)` spans a huge dynamic range at a generic `t` (at `t=1.1`,
  `e^{-11}` vs `e^{-110}` — a ratio of ~10^43), so the probe sees `small ± tiny`
  and the precision ladder climbs to its ceiling without deciding. The
  separation, not the magnitude, is the trigger: rates within ~1 order (−1, −2)
  never provoke it.
- *Profiling nuance (2026-10-16).* A `sample` of the hung process shows the hot
  path is overwhelmingly the **symbolic** engine — `match` / `match_internal_impl`
  / `match_args_internal` (`src/match.c`) plus `apply_replace_repeated_nested` →
  `do_replace_all` (`src/replace.c`) with deeply nested `evaluate_step` recursion
  (`src/eval.c`); MPFR / `zero_test` frames were essentially absent. So in this
  build the stall is dominated by rule-rewrite / pattern-match / candidate churn,
  not the numeric "precision ladder climbs to its ceiling" mechanism described
  above. **Re-diagnose before attempting a fix** — the stated numeric mechanism
  may be stale relative to the current pipeline.

Either way nothing short-circuits the widely-separated-rate case. Note the simp
path has **no dynamic-range guard**, and its deadline machinery
(`g_simp_time_budget` / `g_simp_call_deadline` in `src/simp/simp_search.c`)
defaults to `HUGE_VAL` and is only polled *between* seeds/rounds, never inside a
single transform — so a **bare** `Simplify[...]` search is effectively unbounded
and only a global `TimeConstrained` abort can stop it (hence `$Aborted`, not a
best-so-far result).

**Impact.** Constant-coefficient linear ODE *systems* with a real but spread
spectrum — e.g. `x'=-50x+20y, y'=100x-60y` (eigenvalues −10, −100) — produce a
fundamental-matrix body that is exactly such a sum. **Worked around (M27):**
`dsolve_linsys_tidy` (`src/calculus/dsolve_linsys.c`) routes any
exponential-carrying body through `Expand` rather than `Simplify`; the result is
back-substitution-verified regardless, so the cosmetic loss is harmless. A
direct `Simplify` of a wide-spectrum exponential sum still hangs.

**Suggested direction (core fix, deferred).** Before the numeric ladder, factor
out the dominant exponential (`E^(-10 t) + E^(-100 t) = E^(-100 t)(E^(90 t) + 1)`)
or compare terms in log-magnitude so a term negligible at the sample is not
differenced against a large one; or cap the ladder / return `UNKNOWN` when
per-point dynamic range exceeds the working precision. Given the profiling
nuance above, also investigate the symbolic transform churn (candidate
explosion in the rewrite search), which may be the dominant cost in the current
build regardless of the numeric path.

### Cross-references

- `src/simp/simp_search.c` — the general Simplify pipeline (`simp_search`,
  `g_simp_time_budget`, `g_simp_call_deadline`, `simp_deadline`).
- `src/calculus/dsolve_linsys.c` — `dsolve_linsys_tidy` (the `Expand`-not-
  `Simplify` workaround and its rationale comment).
- `docs/notes/DSOLVE_PLAN.md` — milestone M27 ("Separated-exponent Simplify hang").
- DSolve §2.2.7 corpus systems 636, 650
  (`DSolve_test_status/DE_examples_227.m`).
