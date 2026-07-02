---
description: Produce an implementation plan for a feature before writing any code.
argument-hint: <feature-description>
allowed-tools: Read, Grep, Glob, Bash
---

You are planning a feature. The feature description is `$ARGUMENTS` (or ask
the user if empty). **Do not implement yet.** Your output is a plan.

## Do your homework first

1. Read the project's top-level `README`, `CLAUDE.md` / `AGENTS.md`, and the
   files most relevant to the feature.
2. Run `!git log --oneline -10` to see recent context.
3. Run `!git ls-files '*.md' | head -20` to see what docs exist.

## The plan must cover

### Scope

- What's in scope (bulleted list of the features / changes).
- What's explicitly out of scope.
- Any assumptions you're making about unclear points.

### Files to touch

- List of files you'll create, modify, or delete.
- For each: one sentence on what changes and why.

### Implementation order

A numbered sequence where each step produces a runnable state. Don't plan
a big-bang change; plan increments. Each step should be committable.

### Tests

- Which new tests you'll write.
- Which existing tests you expect to change.
- What the test strategy is (unit, integration, e2e).

### Risks

- What could go wrong during implementation.
- Rollback path if something breaks.

### Estimated size

S / M / L based on file count and behavioral changes. Call out if it
exceeds the standard ≤10-file PR guidance.

## Output format

Use the headers above. Keep it tight — aim for one page. Save the plan so
the user can review it before approving implementation.

## After presenting

Wait for user approval. If they request changes, update the plan. If they
approve, invoke [`/implement-feature`](implement-feature.md) with the plan
as context.

## Related guidelines

- Root [`../AGENTS.md`](../AGENTS.md) engineering principles — change the minimum; extend before create.
- [`../shared/guidelines/pr-and-review.md`](../shared/guidelines/pr-and-review.md) — PR shape expectations.
