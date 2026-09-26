# Fix: `Integrate[Sqrt[t^7/(1-5t^2)] // PowerExpand, t]` Root::conv warning storm + 82 s grind

## Symptom (reported)
`Integrate[t^(7/2)/Sqrt[1-5 t^2], t]` (after PowerExpand) prints an endless stream of
```
Root::conv: refinement did not converge at 298 bits.
Root::conv: refinement did not converge at 431 bits.
Root::conv: refinement did not converge at 564 bits.
```

## Root-cause analysis (verified empirically)
- The integrand is a mixed radical tower; the `Automatic` cascade engages
  `Integrate`ParallelMixedTower` (PMT, `src/internal/mixed/ParallelMixed.m`).
- PMT eventually **declines** (the integral is non-elementary — Mathematica returns a
  `Hypergeometric2F1`; unevaluated is acceptable here) — but it takes **82 s** and prints
  **1560** `Root::conv` lines before returning `Integrate[...]` unevaluated.
- The 1560 warnings are `520` failed `N[Root[...], 80]` calls × 3 precision escalations
  (266 → 399 → 532 bits). They resolve to **two** distinct monic degree-16, all-complex
  (`sturm_real = 0`) minimal polynomials with constant term ≈ 1.55e27, roots `4.4385 ± I/√2`.
- These are re-numericalized **260× each** because PMT's verify-gate samples the candidate
  antiderivative (which contains these Root *constants*) at many points; the Root value is a
  loop-invariant that is recomputed every sample.
- Two genuine defects, both confirmed:
  1. **`root_warn` (`src/root_numeric.c`) writes straight to `stderr`**, bypassing the
     Message subsystem — so the `.m` code's `Quiet[...]` wrappers cannot suppress it and
     there is no repeat-throttle. (Every other emission site uses
     `mth_msg_note_fired(); if (mth_msg_suppressed()) return;` — e.g. `ops_msg` in
     `assoc_ops.c`.)
  2. **`root_numericalize` is not memoised**, so the same deterministic
     `N[Root[p,k], bits]` is recomputed 260×. Verified deterministic: no global MPFR
     state mutation, identical result across 50 isolated repeats.
- NOT in scope (noted follow-up): the MPFR companion-QR solver genuinely fails to refine
  these ill-scaled deg-16 polys where FLINT `qqbar` would succeed. Fixing `N[Root]` to use
  FLINT is a larger, higher-risk change and is unnecessary for correct behaviour here.

## Plan
- [ ] **A. Route `root_warn` through the Message subsystem.** `#include "message.h"`; in
      `root_warn` call `mth_msg_note_fired()` then `if (mth_msg_suppressed()) return;`
      before printing. Matches the established `ops_msg` pattern. → `Quiet[]` now suppresses
      the storm and `Check[]` still sees the firing.
- [ ] **B. Memoise `root_numericalize`.** Rename the MPFR body to
      `root_numericalize_uncached`; add a bounded open-addressed cache (mirror `q2e_cache`)
      keyed by `(expr_hash(root_expr), spec.mode, spec.bits)`, collision-verified by
      `expr_eq`. Cache successes AND failures (NULL). Record whether the uncached call fired
      a message; on a cache hit that recorded a firing, call `mth_msg_note_fired()` (so
      `Check[]` semantics are preserved) but do not reprint. Collapses 1560 solves → 3.
- [ ] Remove temporary `MATHILDA_ROOTDIAG` diagnostics.
- [ ] Verify: repro returns `Integrate[...]` unevaluated, **quietly**, in ~1 s.
- [ ] Regression: root/rootreduce/integrate unit suites; a Charlwood spot-check.
- [ ] Bump `src/version.h` 0.204 → 0.205; changelog note in
      `docs/spec/changelog/2026-09-21.md` (Monday of this ISO week).

## Review — DONE (v0.205)
- **A** `src/root_numeric.c` `root_warn`: added `#include "message.h"`; now
  `mth_msg_note_fired(); if (mth_msg_suppressed()) return;` before printing.
- **B** `src/root_numeric.c`: renamed the MPFR entry to `root_numericalize_uncached`; added a
  256-slot open-addressed memo (`g_rnum_cache`) keyed by `(expr_hash(Root), mode, bits)`,
  collision-verified by `expr_eq`, caching successes and failures, recording whether the
  uncached call fired a message and re-noting it on hits (preserves `Check[]`).
- **C (found mid-fix)** `src/linalg/inv.c`: the same integral produced a SECOND storm,
  `Inverse::matsq: Argument {} …` (574 lines, present all along — PMT calls `Inverse[{}]` in a
  loop inside `Quiet[]`), because inv.c's 5 `Inverse::sing`/`::matsq` sites also wrote bare to
  `stderr`. Added `#include "message.h"` + a static `inv_warn` (note_fired + suppressed) and
  routed all 5 through it.
- Removed the temporary `MATHILDA_ROOTDIAG` diagnostics.
- **Verified:** repro → `Integrate[t^(7/2)/Sqrt[1-5 t^2], t]` unevaluated, **0** warnings of
  either kind, 82 s → 27 s. Profiling confirms 0 samples now in `root_numericalize`; the
  residual is deg-16 `ToNumberField`/`qqbar` (pre-existing PMT gap, out of scope).
- **Semantics:** unquieted `N[Root[bad-k]]` still prints; `Quiet[…]` suppresses; `Check[…]`
  returns the fail-expr. `N[Root]` values unchanged (√2, plastic number, deg-16 pair).
- **Regression:** `root_numeric`/`rootreduce`/`nroots`/`findroot`/`integrals`/
  `parallelmixedtower`/`linalg`/`linearsolve` suites pass; Charlwood P2/P4/P8/A2/A3/P9 verify
  (residual 0). `make check-c99` exit 0; no build warnings.
- **Not done (noted):** (1) MPFR companion-QR can't refine ill-scaled deg-16 polys where FLINT
  `qqbar` would — a larger, higher-risk `N[Root]` robustness change; unnecessary for correct
  behaviour here. (2) ~70 other bare `fprintf(stderr,"Head::tag:…")` sites across `src/` share
  the Quiet-bypass class (Det::matsq, MatrixPower::matsq, Dot::dotsh, …) — a systemic hygiene
  sweep, not done here. (3) The escalation loop reprints an out-of-range-`k` message 3× within
  one call (pre-existing, benign, now memoised across calls).
- **Pending:** version bumped to 0.205 + changelog written; git commit/tag left to the user.
