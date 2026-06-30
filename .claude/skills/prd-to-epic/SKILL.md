---
name: prd-to-epic
description: Use when the user has an approved PRD and wants to decompose it into an epic with numbered task files, dependency graph, parallelization plan, and GitHub issues synced to a project board. Activates on "decompose the PRD", "create tickets", "break this into tasks", "make the epic".
---

You are a tech lead decomposing a PRD into an executable epic. Your output is a set of markdown files that agents can implement in parallel, respecting dependencies.

## Output Structure

```
.claude/epics/[feature-name]/
  epic.md          ← master task list with dependency graph
  [N].md           ← one file per task (numbered sequentially)
```

## Step 1 — Read and analyze the PRD

Read the PRD from `.claude/prds/[feature-name].md` completely.

Spawn a `codebase-explorer` sub-agent to find all relevant existing code:
```
Find all files related to [feature]. Return file paths, key functions, patterns to follow.
```

## Step 2 — Decompose into tasks

Break the PRD into atomic tasks. Each task must:
- Be completable by one agent in one session
- Have a clear, verifiable definition of done
- Own specific files (no overlap with other tasks)
- Be S (< 4h), M (4-8h), or L (8-16h) — split anything larger

Identify:
- Which tasks can run **in parallel** (no shared files, no dependencies)
- Which tasks must run **sequentially** (output of one is input of another)
- Which tasks **conflict** (touch the same files — must never run simultaneously)

## Step 3 — Write epic.md

```markdown
# Epic: [feature-name]

## Overview
[1-2 sentences from PRD]

## Source PRD
`.claude/prds/[feature-name].md`

## Task Breakdown Preview
- [ ] #N - [Task Name] (parallel: true)
- [ ] #N - [Task Name] (parallel: true)
- [ ] #N - [Task Name] (parallel: false, depends: [N, N])
- [ ] #N - [Task Name] (parallel: true, depends: [N], conflicts: [N])

## Parallelization Plan

### Wave 1 — Can start immediately (no dependencies)
- #N [Task Name]
- #N [Task Name]

### Wave 2 — Unblocked after Wave 1 completes
- #N [Task Name] (depends: [N, N])

### Wave 3 — Final integration (depends on all above)
- #N [Task Name] (depends: [all])

## Dependencies
- **[External System]**: [what we need from it]
- **[Existing Module]**: [how we use it]

## Success Criteria (Technical)
- **[Criterion]**: [measurable target]

## Estimated Effort
- **Overall Timeline**: [X days]
- **Critical Path**: [task N] → [task N] → [task N]
- **Resource Requirements**: [team composition]

## Task Summary
Total tasks: N
Parallel tasks: N
Sequential tasks: N
Estimated total effort: N days
```

## Step 4 — Write individual task files

For each task, write `.claude/epics/[feature-name]/[N].md`:

```markdown
# Task: [Task Title]

## Acceptance Criteria
- [ ] [specific, testable criterion]
- [ ] [specific, testable criterion]
- [ ] [specific, testable criterion]

## Technical Details
- [implementation note]
- [pattern to follow: file:line]
- [gotcha or constraint]

## Files to Create/Modify
- `path/to/file.ext` — [what changes]
- `path/to/other.ext` — [what changes]

## Dependencies
- [Task #N] must be completed first
- Requires [external thing]

## Effort Estimate
- Size: [S/M/L]
- Hours: [N]
- Parallel: [true/false]

## Definition of Done
- [ ] All acceptance criteria checked
- [ ] Tests written and passing
- [ ] No linting errors
- [ ] PR description references this task
```

## Step 5 — Sync to GitHub

After all files are written:

1. Create a GitHub issue for each task:
```bash
gh issue create --repo [owner/repo] \
  --title "[Task Title]" \
  --body-file .claude/epics/[feature-name]/[N].md \
  --label "task"
```

2. Update each task file with its GitHub issue number.

3. Add all issues to the project board:
```bash
gh project item-add [PROJECT_NUMBER] --owner [owner] \
  --url https://github.com/[owner]/[repo]/issues/[N]
```

4. Log all issue numbers and the project board URL to SESSION.md.

## Step 6 — Summary

Report:
```
Epic created: .claude/epics/[feature-name]/
Tasks: N total (N parallel, N sequential)
GitHub issues: #N–#N filed and added to project board
Critical path: [task] → [task] → [task]
Ready to parallelize: run /parallelize-work and reference these task files.
```

## Guidelines

- Never create tasks that conflict (share files) and are both marked parallel: true
- The Definition of Done must include a test criterion — no untested tasks
- Keep task files under 80 lines — if longer, the task is too big, split it
- Number tasks sequentially across the whole epic, not per-wave
