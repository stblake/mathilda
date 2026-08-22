---
ticket: DEMO-2
created: 2026-08-21T00:00:00-04:00
researcher: Michael Sollami
topic: "NMinimize returns constrained solutions that violate their constraints by ~1e-4"
type: research
lifecycle: active
full_research: thoughts/shared/tickets/DEMO-2/research.md
---

# Research Summary: NMinimize constrained-feasibility bug

**Full research (appendix)**: `thoughts/shared/tickets/DEMO-2/research.md`

## Recommendation

Keep `nm_eval_pen` returning the squared penalty. Define the feasibility tolerance on
the *actual* violation and square it at the definition site, so `NM_FEAS_EPS` and
`NM_FEAS_FINAL` mean what their comments claim without changing the quantity six
engines do arithmetic on. Add feasibility assertions that are two-sided, since every
existing one is a one-sided band a violating point already passes.

## Options Considered

1. **Return a violation norm from `nm_eval_pen`** — dimensionally clean, and the first
   choice. Refuted by research: 6 of 8 engines compute `f + 1e6*penalty` and feed it to
   Metropolis acceptance functions, NelderMead's simplex ranking, and DIRECT's convex
   hull. Rescaling the penalty changes what those five algorithms search.
2. **Square the tolerance at the definition site** — arithmetic engines bit-identical,
   relationship explicit in code, cannot drift. **Chosen.**
3. **`sqrt(p)` at the four threshold sites only** — equally safe for the engines, but
   repeats the sqrt at every site and puts one in DE's hot loop.
4. **Proceed with the norm and re-tune the five engines** — turns a bug fix into an
   optimizer rewrite and invalidates benchmarks 79-86.

## Decision Criteria

- **Blast radius.** `nm_eval_pen` has exactly one caller, so its consumers are fully
  enumerable — and enumerating them is what disqualified option 1.
- **The squared scale is load-bearing.** `NM_PENALTY_MU = 1e6` is calibrated against it.
- **The defect is readability, not representation.** The constant's units were implicit;
  making them explicit fixes the cause without disturbing a convention six sites rely on.
- **Speed is checkpointed, not assumed.** Correctness wins, but a >~3x regression on the
  constrained cases stops the work for a re-decision.

## Open Questions

### Unresolved

_None._

### Resolved

- [x] Which fix? — Square the tolerance at the definition site (revised after research
      refuted the norm approach).
- [x] Fix `NM_FEAS_FINAL` too? — Yes, same treatment, violation tolerance chosen
      deliberately rather than inherited.
- [x] Acceptable speed cost? — Measure and report; stop if worse than ~3x.
- [x] Risk of flipping to `Infinity`? — Keep the final threshold loose; six named tests
      must be re-run.
- [x] Is a feasibility test really missing? — Yes, structurally: every existing check is
      one-sided.

## Requires Approval

None outstanding. The fix approach was re-confirmed with the human after research
contradicted the initial choice, and the speed checkpoint is already agreed.
