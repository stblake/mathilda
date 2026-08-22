---
ticket: DEMO-2
created: 2026-08-21
type: plan
lifecycle: active
status: draft
full_plan: thoughts/shared/tickets/DEMO-2/plan.md
---

# Plan Summary: NMinimize constrained-feasibility fix

**Full plan (appendix)**: `thoughts/shared/tickets/DEMO-2/plan.md`

## Recommendation

State each feasibility tolerance on the *actual* constraint violation and square it at the
definition site in `findmin_internal.h`. Leave `nm_eval_pen`'s returned quantity untouched,
so the six engines that compute `f + NM_PENALTY_MU * p` are bit-identical. Add two-sided
feasibility assertions to `tests/test_nminimize.c`, which currently has none.

## Options Considered

1. **Return `sqrt(total)` from `nm_eval_pen`** — dimensionally cleanest, and the approach
   originally chosen. Refuted during research: 6 of 8 engines feed the penalty into
   Metropolis acceptance, simplex ranking, or convex-hull slopes, so rescaling it changes
   what they search.
2. **Square the tolerance at the definition site** — engines untouched, units explicit in
   source, cannot drift. **Chosen.**
3. **`sqrt(p)` at the four threshold sites** — equally safe, but four sites to keep in sync
   and a `sqrt` in DE's inner loop.
4. **Norm + re-tune the five affected engines** — turns a bug fix into an optimizer rewrite
   and invalidates benchmarks 79–86.

## Decision Criteria

- `nm_eval_pen` has exactly one caller, so its consumers are fully enumerable — and
  enumerating them is what disqualified option 1.
- `NM_PENALTY_MU = 1e6` is calibrated against the squared scale; that scale is load-bearing.
- The defect is *implicit units*, not the representation. Make the units explicit.
- Speed is checkpointed at ~3x, not assumed.

## Decisions

- Square the tolerance at the definition site rather than normalise the penalty.
- `NM_FEAS_FINAL` keeps its current numeric value (`1e-6` = `(1e-3)^2`) so the
  Infinity/give-up path provably cannot regress.
- Correctness over speed, with a hard stop-and-report at ~3x regression.

## Non-goals

- Not changing `nm_eval_pen`'s return value or re-tuning any engine.
- Not adding an infeasibility warning (follow-up, recorded in the bug note).
- Not fixing the `PenaltyFunction` tolerance-semantics wart.
- Not touching `FindMinimum`'s local solver — research proved it independent.

## Open Questions

### Unresolved

_None._

### Resolved

- [x] Fix approach — square at the definition site (revised after research refuted the norm).
- [x] `NM_FEAS_FINAL` too — yes, effective value held constant.
- [x] Speed cost — measure and report, stop at ~3x.
- [x] Fallback if the checkpoint trips — decouple DE's convergence gate (`nm_de.c:159`).
- [x] Is a feasibility test missing — yes, every existing check is one-sided.

## Requires Approval

_None._ Approach, checkpoint, and fallback agreed before the plan was written.

## Architecture Impact

- New services introduced: none
- APIs changed: none
- Data crossing a service boundary: none
- New external dependency: none
- Deployment topology change: none

## Subsystems & Dependencies

- Subsystems touched: none declared
- Interdependencies surfaced: none
