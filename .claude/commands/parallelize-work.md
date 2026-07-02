---
description: Decompose a feature or epic into parallel and sequential work streams, assign file ownership, and optionally dispatch the parallel-worker agent to execute them.
argument-hint: <feature-description-or-plan-path>
allowed-tools: Read, Grep, Glob, Bash
---

Analyze `$ARGUMENTS` (a feature description or path to a plan file) and
decompose it into parallel and sequential work streams with file ownership
assignments. Optionally execute via the `parallel-worker` agent.

## Step 1 — Understand the work

If `$ARGUMENTS` is a file path, read it. Otherwise treat it as a description
and ask for clarification on scope if needed.

Read the relevant parts of the codebase to understand:
- Which files will be touched
- Which changes have hard dependencies on others
- Which changes are truly independent

## Step 2 — Identify streams

Classify each sub-task:

| Stream | Work | Files | Parallel? | Depends on |
| --- | --- | --- | --- | --- |
| A | <description> | <file list> | yes | — |
| B | <description> | <file list> | yes | — |
| C | <description> | <file list> | no | A, B |

Rules:
- **Parallel** = no shared files, no logical dependency.
- **Sequential** = must wait for another stream's output (schema, type, API contract).
- A stream that only reads a file owned by another stream is still parallel —
  reads don't conflict.

## Step 3 — Check for conflicts

Verify no two parallel streams touch the same file. If they do, either:
1. Split the file's changes and assign different line ranges, or
2. Make one stream sequential (wait for the other).

## Step 4 — Write the stream specs

For each stream, produce a spec block:

```
## Stream A — <name>
Branch: feat/<scope>-stream-a
Files: <exact list>
Work:
  - <acceptance criterion 1>
  - <acceptance criterion 2>
Parallel: yes
Depends on: —
Estimated size: S / M / L
```

## Step 5 — Present for approval

Show the decomposition. Ask:
- Does the stream breakdown look right?
- Are the file assignments complete?
- Should any streams be merged or split?

## Step 6 — Execute (if approved)

If the user approves, either:

**Option A — Dispatch parallel-worker agent** (recommended for 3+ streams):
Hand the decomposition to the `parallel-worker` agent with the stream specs.
It creates worktrees, spawns sub-agents, and consolidates results.

**Option B — Guide sequential execution** (for 1–2 streams or when the user
wants visibility into each step):
Walk through each stream in order, running tests after each one.

## Output format

```
## Work Decomposition: <feature name>

### Summary
- Total streams: N
- Can run in parallel: N
- Sequential (blocked): N
- Estimated total: <S/M/L>

### Stream breakdown
[table from step 2]

### Execution plan
Week 1: Streams A, B (parallel)
Week 2: Stream C (after A and B merge)

### File ownership map
[who owns what]

### Risks
- <anything that could cause the decomposition to fail>
```

## Related

- Agent: [`../agents/parallel-worker.md`](../agents/parallel-worker.md)
- Guideline: [`../shared/guidelines/worktree-workflow.md`](../shared/guidelines/worktree-workflow.md)
- Command: [`plan-feature.md`](plan-feature.md) — use first for single-stream work.
