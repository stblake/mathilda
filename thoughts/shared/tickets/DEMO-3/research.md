---
ticket: DEMO-3
created: 2026-08-22T00:00:00-04:00
researcher: Michael Sollami
source_sha: 70fb7bdd985279afceb1275c97124eb3b4ac8f4f
branch: claude/where-are-you-5f1147
repository: mathilda
topic: "Audit of tolerance-hiding in the NMinimize test suite"
tags: [research, audit, tests, tolerances, nminimize, DEMO-3]
subsystems: [numerical_calculus]
type: research
lifecycle: active
status: complete
last_updated: 2026-08-22
last_updated_by: Michael Sollami
---

# Research: Tolerance-hiding in the NMinimize test suite

**Date**: 2026-08-22
**Git Commit**: 70fb7bdd985279afceb1275c97124eb3b4ac8f4f
**Branch**: claude/where-are-you-5f1147
**Ticket**: DEMO-3

## TL;DR

Two known instances were not a coincidence. Of 83 tests, **32 constrained tests never
check their constraints at all**, and mutation proves **57 accept a point 10% wrong**
while only 23 miss a 10%-wrong objective. The suite tests the objective, not the answer.
Nothing was tuned to green after the fact — the looseness was authored in, under a stated
policy that covers objectives and is silent on feasibility.

## Summary

Two independent methods agree.

**Structural.** Classifying all 83 tests by assertion shape gives 32 constrained tests
with no feasibility assertion of any kind, 2 with feasibility but no objective, 17
one-sided feasibility ceilings, and 5 shape-only tests.

**Mutation.** Perturbing `nm_build_result` to return a deliberately wrong answer, with
the suite instrumented to continue past failures, measures how wrong an answer each test
will accept. A point 10% wrong is caught by 26 of 83; an objective 10% wrong by 60 of 83.
That asymmetry *is* the finding: the suite is roughly twice as good at noticing a wrong
number as a wrong answer.

The two agree on mechanism: every test that asserts only the objective is blind to the
returned point, and every test that asserts the point or a constraint catches it. So the
original bug was not bad luck. It was the predictable consequence of a suite whose stated
tolerance policy governs objectives and says nothing about feasibility.

Git history exonerates the authors of the thing one would suspect: across 37 commits
there is **no** instance of a tolerance being widened to make a failing test pass.

## Open Questions

### Unresolved

_None._

### Resolved

- [x] **What counts as a finding?** — Structural classification of all 83 by assertion
  shape, plus mutation evidence on the suspects. Both were produced.
- [x] **Is it in scope to change what tests assert?** — Yes: tighten numbers where the
  problem supports it, and add an assertion only where shape analysis proves a test is
  blind to a whole class.
- [x] **How to handle the integer gap?** — Extended the instrument with
  `MATHILDA_MUTATE_INT`, which flips one integer coordinate. Now measured, not assumed.
- [x] **Were tolerances tuned to green?** — No. 37 commits, essentially all additions;
  the only tolerance *removals* in the whole history are from this branch's own
  `d4eb10fa`. Looseness was authored in, not retrofitted.

## Requires Approval

_None._ Scope, finding bar, and fix policy were all settled before this document.

---

## Research Question

> Ticket DEMO-3. Audit the NMinimize test suite (tests/test_nminimize.c) for
> tolerance-hiding: tests whose tolerances are loose enough to pass against a wrong
> answer. Two instances are known. I need to know how widespread the pattern is, then fix
> what's hiding failures.

## Detailed Findings

### 1. Mutation results (the hard evidence)

Instrument: `nm_build_result` (`findmin_nm_common.c`) gated on three env vars —
`MATHILDA_MUTATE_PT` (perturb the returned point, leave the objective true),
`MATHILDA_MUTATE_OBJ` (the converse), `MATHILDA_MUTATE_INT` (flip one integer
coordinate). `ASSERT_STR_EQ` was temporarily made non-fatal so one run yields the whole
pass/fail map instead of stopping at the first failure. Baseline: 0 failures, 83 tests.

**Point wrong, objective right** — the shape of the original bug:

| relative error in the point | tests that caught it | blind |
|---:|---:|---:|
| 1e-6 | 4 | 79 |
| 1e-4 | 10 | 73 |
| 1e-2 | 24 | 59 |
| **1e-1** | **26** | **57** |

**Objective wrong, point right** — the converse:

| relative error in the objective | caught | blind |
|---:|---:|---:|
| 1e-4 | 23 | 60 |
| 1e-2 | 51 | 32 |
| **1e-1** | **60** | **23** |

**One integer coordinate flipped by ±1:**

| mutation | caught |
|---|---|
| +1 | 7 — `3sat_feasibility`, `adjacency_assignment`, `cardinality_portfolio`, `max_independent_set`, `multiknapsack`, `qap_assignment`, `sudoku_latin` |
| −1 | 6 — as above minus `max_independent_set`, plus `txncost_portfolio` |

Integer-heavy tests blind to a flip: `integer_domain_value`, `integer_domain_heads`,
`integer_domain_list`, `integer_domain_alternatives`, `mixed_integer`,
`job_scheduling`, and `txncost_portfolio` (to +1 only).

**Instrument caveat, stated because it nearly produced a false headline.** The first INT
sweep reported 0 of 83 — which would have been a spectacular finding and was an artefact
of a stale CMake object. The standalone binary applied the mutation correctly while the
test binary did not. Re-running after a forced rebuild gave 7. The PT and OBJ sweeps were
re-run on the known-good build and reproduced their earlier numbers exactly.

### 2. Structural classification (all 83 tests)

**List A — constrained, NO feasibility assertion at all (32 tests).** The headline class.
Grouped by why they exist:

- *Early constraint suite, objective only:* `disk_linear`, `quadratic_linear`,
  `linear_program`, `equality_constraint`, `chained_inequality`, `equation_system`
- *Integer/domain tests:* `integer_domain_value`, `integer_domain_heads` (no objective
  either — shape only), `mixed_integer`, `integer_domain_alternatives`,
  `integer_domain_list`, `indexed_table_constraints`, `indexed_rosenbrock`,
  `indexed_real_coefficient`, `region_expansion_rescue`
- *Method/sub-option regressions on box-constrained problems:* `sa_suboptions`,
  `griewank_simulatedannealing`, `sa_deceptive_landscapes`,
  `schaffer2_simulatedannealing`, `griewank_differentialevolution`,
  `griewank_neldermead`, `randomsearch_searchpoints_verbatim`, `bukin6_no_warning`,
  `initial_points`, `penalty_function`, `autocompile_parity_and_fallback`,
  `symbol_indirection`, `search_points_honored`, `de_options_effective`,
  `de_boundary_no_stagnation`
- *NMaximize wrappers:* `nmaximize_constrained`, `min_max_duality`
- *Applied testbed:* `gaussian_well`, `modified_ackley`, `job_scheduling`

**List B — feasibility but no objective (2).** `returned_point_feasible`,
`returned_point_feasible_all_methods`. Converse blindness: they confirm the point is
feasible without ever checking it is optimal, so a feasible-but-useless answer passes.

**List C — one-sided feasibility ceilings (17).** Slacks: `9.001` (1e-3), `25.001`
(1e-3), `3.0001` (1e-4), `100.0001` (1e-4), `5.1201` (1e-4), `100.001`/`200.001`
(1e-3), `15.001` (1e-3), `0.027001` (1e-6), `1.0001` (1e-4), `5.001` (1e-3), plus
several `>= -1e-6` non-negativity floors. Each catches overshoot and is blind to a point
that never *reaches* the boundary — the direction the original bug went.

**List E — shape-only (5).** `result_shape`, `options_nonempty`,
`max_iterations_accepted`, `arity_error_unevaluated`, `integer_domain_heads`.
Legitimately non-numeric except the last, which is a constrained solve asserting only that
the result has head `Integer`.

### 3. Principled vs. tuned-to-green

The file's own top-of-file docstring states a policy: *"Objective tolerances are looser
than FindMinimum's (global heuristics + penalty constraints), tight where the problem is
convex/linear."* That policy is defensible and observably followed —
`optimal_liquidation` and `minimax_chebyshev` (convex) assert at 1e-4/1e-3, while
stochastic global searches sit at 1e-2.

**The policy has no feasibility clause.** That is the whole defect. Feasibility is not an
approximation subject to search variance — a point either satisfies its constraints or it
does not — but the suite treats it with the same discretion it applies to objectives, and
in 32 cases with no discretion at all because the check is simply absent.

Git history confirms this is an authoring gap rather than erosion: 37 commits touching the
file, essentially all additive, no tolerance ever widened to rescue a failing test.

## Results after the fix (measured, not predicted)

Re-run on a forced-clean build after Phase 2/3. The plan originally predicted >= 34 and
>= 36; the reviewer flagged those as invented and they were demoted to "report what is
measured" before implementation. **Both would have been missed** — which is the point of
having demoted them.

| probe | before | after |
|---|---:|---:|
| PT -1e-4 | 10 | **16** |
| PT -1e-2 | 24 | **30** |
| PT -1e-1 | 26 | **32** |
| INT +1 | 7 | **8** |
| unmutated | 0 failures | 0 failures |

INT moved by only one because the new integer-domain assertions are `x + 2y >= 3` style
lower bounds, and flipping a coordinate *upward* keeps them satisfied. A one-directional
constraint is one-directionally sensitive — the same lesson as the one-sided ceilings.

Structural state after the change: every constrained test now either asserts its
constraints or carries a comment saying why it does not. 18 deliberate omissions are
documented (16 box-only, 1 unconstrained-and-misfiled, 1 deferred with reason).

## Code References

- `tests/test_nminimize.c:1-24` — the stated tolerance policy
- `tests/test_nminimize.c:106-133` — the six objective-only constrained tests
- `tests/test_nminimize.c:140-143` — feasibility without objective
- `tests/test_nminimize.c:220-249` — integer-domain tests, objective + head only
- `src/numerical_calculus/findmin_nm_common.c` `nm_build_result` — mutation site
- `tests/test_utils.h` `ASSERT_STR_EQ` — made non-fatal for the sweep

## Architecture Insights

A test suite has two independent things to assert about an optimiser's answer: **is this
the right value** and **is this a legal point**. They fail independently, they need
different tolerances, and only the first is a matter of degree. This suite has a
thought-through policy for the first and no vocabulary for the second — so feasibility got
handled ad hoc, per test, by whoever wrote it, which produced 32 omissions and 17
one-sided ceilings without anyone making a bad decision.

The mutation asymmetry (60 catch a wrong objective, 26 a wrong point) is that absence,
measured.

## Historical Context (from thoughts/)

- `NMINIMIZE_FEASIBILITY_BUG.md` — the original bug and its "second instance" section,
  which this audit was commissioned from.
- `thoughts/shared/tickets/DEMO-2/research.md` — established that every existing
  feasibility check is one-sided.

## Related Research

- `thoughts/shared/tickets/DEMO-2/plan.md` — the two-threshold fix whose test gap
  prompted this audit.
