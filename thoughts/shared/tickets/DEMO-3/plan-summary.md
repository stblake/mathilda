---
ticket: DEMO-3
created: 2026-08-22
type: plan
lifecycle: active
status: draft
full_plan: thoughts/shared/tickets/DEMO-3/plan.md
---

# Plan Summary: NMinimize test-suite tolerance audit fixes

**Full plan (appendix)**: `thoughts/shared/tickets/DEMO-3/plan.md`

## Recommendation

Give feasibility its own stated policy in the file header, add feasibility assertions to
the tests with real constraints, and document as deliberate the omissions that are
vacuous. Measure the result with the same mutation instrument that found the problem, then
remove the instrument.

## Options Considered

1. **Tighten numbers only** — cannot fix the class; a test that never asserted feasibility
   stays blind however tight its objective bound gets.
2. **Policy + assertions where constraints are real + documented omissions** — **chosen.**
3. **Uniform assertion on all 32** — half are box-only where the solver clamps, so the
   assertion cannot fail; a check that cannot fail is the same theatre being audited.

## Decision Criteria

- The defect is a missing category, not carelessness — so state the category.
- A check that cannot fail is worse than an acknowledged gap.
- Assert against measured behaviour with margin, never at a threshold constant.

## Decisions

- Feasibility policy written into the file header.
- Real constraints get assertions; box-only tests get an explicit comment.
- Re-measure with the audit's own instrument; remove it last.

## Non-goals

- No solver changes; the instrument is removed before commit.
- Objective tolerances that are loose for a stated reason stay loose.
- `test_fixed_charge_flow`'s DEMO-2 regression is a separate ticket.
- Other test files not audited.

## Open Questions

### Unresolved

_None._

### Resolved

- [x] Finding bar — structural classification of all 83 plus mutation.
- [x] Fix scope — assertions may be added where shape analysis proves blindness.
- [x] Integer measurement gap — closed with `MATHILDA_MUTATE_INT`.
- [x] Box-only treatment — comment, not a vacuous assertion.

## Requires Approval

_None._

## Architecture Impact

- New services introduced: none
- APIs changed: none
- Data crossing a service boundary: none
- New external dependency: none
- Deployment topology change: none

## Subsystems & Dependencies

- Subsystems touched: none declared
- Interdependencies surfaced: none
