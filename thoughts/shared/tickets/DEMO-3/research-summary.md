---
ticket: DEMO-3
created: 2026-08-22T00:00:00-04:00
researcher: Michael Sollami
topic: "Audit of tolerance-hiding in the NMinimize test suite"
type: research
lifecycle: active
full_research: thoughts/shared/tickets/DEMO-3/research.md
---

# Research Summary: Tolerance-hiding in the NMinimize test suite

**Full research (appendix)**: `thoughts/shared/tickets/DEMO-3/research.md`

## Recommendation

The pattern is systemic, not two incidents. Fix it by giving feasibility its own policy
rather than by adjusting numbers: add a feasibility assertion to the constrained tests that
have none, make one-sided ceilings two-sided where the constraint is an equality, and where
a bound genuinely cannot be tightened, say so in the test with the measured reason instead
of leaving a loose number that reads as an assertion.

## The finding, in one table

| probe | tests that catch it | blind |
|---|---:|---:|
| returned point 10% wrong, objective right | 26 / 83 | **57** |
| objective 10% wrong, point right | 60 / 83 | 23 |
| one integer coordinate flipped ±1 | 7 | 76 |

| structural class | count |
|---|---:|
| constrained, **no feasibility assertion at all** | **32** |
| feasibility asserted but objective never checked | 2 |
| one-sided feasibility ceilings | 17 |
| shape-only, no numeric content | 5 |

The suite is about twice as good at noticing a wrong *number* as a wrong *answer*.

## Options Considered

1. **Tighten numbers only** — small, reviewable, and cannot fix the actual class: a test
   that never asserted feasibility still won't, however tight its objective bound gets.
2. **Add feasibility assertions where shape analysis proves blindness, and tighten where
   the problem supports it** — targeted by the audit's own evidence. **Chosen.**
3. **Blanket-tighten everything** — would break the many tests whose looseness is
   principled (stochastic global search genuinely varies), and would encode today's
   measurements as contracts.

## Decision Criteria

- **Two methods agreed.** Structural classification and mutation independently identify
  the same tests, so the finding does not rest on either alone.
- **The defect is a missing category, not sloppiness.** The file states a tolerance policy;
  it covers objectives and is silent on feasibility. 37 commits, no tolerance ever widened
  to rescue a failing test.
- **Feasibility is not a matter of degree.** An objective may legitimately be loose because
  the search is stochastic. A point either satisfies its constraints or it does not.

## Open Questions

### Unresolved

_None._

### Resolved

- [x] What counts as a finding? — Structural classification of all 83, plus mutation on the
      suspects. Both produced.
- [x] May I change what tests assert? — Yes, where shape analysis proves blindness.
- [x] The integer gap? — Closed by extending the instrument to flip an integer coordinate.
- [x] Were tolerances tuned to green? — No. Authored loose, never widened afterwards.

## Requires Approval

_None._ Finding bar, fix scope, and instrument coverage were all agreed before this
document was written.
