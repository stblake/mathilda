---
description: Execute an approved feature plan end-to-end through PR.
argument-hint: <feature-name-or-plan-path>
allowed-tools: Bash, Read, Edit, Write, Grep, Glob
---

You are implementing an approved plan. If `$ARGUMENTS` names an existing
plan file, read it. Otherwise ask the user for the plan or invoke
[`/plan-feature`](plan-feature.md) first — **do not improvise an
implementation without a plan.**

## Execute the plan in order

Follow the plan's numbered steps. After each step:

- Run the relevant tests.
- Commit if the step produces a stable state and the project permits
  multi-commit PRs.
- If you discover the plan is wrong, stop and update the plan with the
  user before continuing.

## Implementation rules

- **Extend, don't create.** Add new files only when the plan says to, or
  when there's no natural home for the change.
- **Match existing patterns.** Don't introduce new state libraries, test
  runners, or style approaches mid-feature.
- **Minimum change.** No unrelated cleanup. No speculative abstraction.
- **Tests live with the change.** If you add behavior, you add a test. If
  you change behavior, you update the test.

## Verify before you declare done

- [ ] All tests pass.
- [ ] Typechecker / linter clean.
- [ ] Feature works against a sample request (UI: drive it in the browser).
- [ ] No committed secrets, no committed console/print statements.
- [ ] Commit history tells a story — squash WIP noise.

## Open the PR

Invoke [`/review-pr`](review-pr.md) against your own branch first as a
self-review. Then open the PR with:

- Title: imperative, ≤ 70 chars.
- Summary: 1–3 bullets.
- Test plan: checklist a reviewer can run.
- Out of scope: things noticed but not fixed.

See [`../shared/guidelines/pr-and-review.md`](../shared/guidelines/pr-and-review.md).

## If you get stuck

Don't push through silently. Report the blocker to the user with:

- What step of the plan you're on.
- What you tried.
- What the failure mode is.
- Two or three ways forward with tradeoffs.

## Related assets

- Skill: [`../skills/python-service-implementation/SKILL.md`](../skills/python-service-implementation/SKILL.md)
- Skill: [`../skills/react-feature-delivery/SKILL.md`](../skills/react-feature-delivery/SKILL.md)
- Skill: [`../skills/terraform-module-workflow/SKILL.md`](../skills/terraform-module-workflow/SKILL.md)
