---
description: Work a production issue from observation to reproduction to root cause to fix.
argument-hint: <symptom-or-ticket>
allowed-tools: Bash, Read, Grep, Glob, WebFetch
---

A production issue is reported. Symptom / context is `$ARGUMENTS`. Your job
is to find and fix the root cause — not just make the symptom go away.

## Phase 1 — Observe

1. What exactly is broken? Stack trace, error rate, user-visible behavior,
   time it started.
2. What changed recently? `!git log --oneline --since='2 days ago'`
3. Is it a full outage, degraded, or intermittent?
4. What's the blast radius? Which users, regions, or features are affected?

If any of these are unknown, ask the user before going further.

## Phase 2 — Reproduce

- Can you trigger it locally? Same inputs, same config?
- If not local, is there a staging or shadow environment?
- If the bug is intermittent, what state makes it more likely? Time of day,
  load, specific inputs?

**Don't theorize past this phase** without a reproduction. A fix without
a repro is a guess.

## Phase 3 — Isolate

- Bisect recent commits if the bug is new: `!git bisect start`.
- Narrow the input space to the minimum that triggers the failure.
- Identify the *layer* responsible — is it the client, the service, the
  database, the cloud?

## Phase 4 — Root cause

Answer *why* the code was wrong, not just *what* to change:

- Was there a missing check? Why wasn't it there?
- Was there a race? What ordering assumption broke?
- Was there a config / env difference between prod and elsewhere?
- Was there a dependency update?

## Phase 5 — Fix

- Write a failing test that reproduces the issue.
- Make it pass with the minimum change.
- Verify the full test suite still passes.
- See [`../skills/bug-triage-and-debugging/SKILL.md`](../skills/bug-triage-and-debugging/SKILL.md).

## Phase 6 — Land safely

- Open a PR labeled as a hotfix if the project distinguishes.
- PR description explains: symptom, root cause, fix, regression test.
- If the bug is severe, coordinate with the on-call / incident channel —
  don't just land the fix and leave.
- Follow-ups: did we learn something about monitoring? File a ticket.

## Postmortem contribution

For anything user-visible or that caused an outage, the fix PR is the
start — not the end. Capture:

- Timeline (when it started, when it was noticed, when it was fixed).
- Impact (scope, duration).
- Root cause (one clear sentence).
- What made it hard to catch.
- Follow-up actions.

Don't write the postmortem in this library — write it in the project.
This command just reminds you to.

## Guardrails

- Don't `--no-verify` or `.skip` tests to land the fix.
- Don't land an untested fix to an unreproducible bug.
- Don't revert a large commit when the actual bug is one line of it —
  narrow the revert.
