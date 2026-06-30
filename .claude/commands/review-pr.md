---
description: Structured review of a PR or branch diff — correctness, security, design, readability.
argument-hint: [PR# or branch]
allowed-tools: Bash, Read, Grep, Glob
---

Review the code under review. If `$ARGUMENTS` is a PR number, fetch it via
`gh pr view` and `gh pr diff`. If it's a branch name, diff against the
project's default branch. If empty, review the current branch against its
upstream.

## Gather context

!`gh pr view $ARGUMENTS 2>/dev/null || git log --oneline origin/main..HEAD`
!`gh pr diff $ARGUMENTS 2>/dev/null || git diff origin/main...HEAD --stat`

## Review buckets

Work through each bucket; report findings grouped by **blocker / warning / nit**.

### Correctness

- Does it do what the PR says? Read the description and confirm the diff
  matches.
- Edge cases: empty inputs, error paths, concurrency if relevant.
- Tests assert on behavior, not implementation? Would they catch a real regression?

### Security

Check against [`../shared/guidelines/security-baseline.md`](../shared/guidelines/security-baseline.md):

- Secrets in source or logs.
- Unescaped user input into SQL / HTML / shell.
- New public endpoints without auth or rate limiting.
- Overly broad IAM / security-group rules.
- Dependencies: freshly published, obscure, or transitively pulling in heavy libs.

### Design

- Consistent with codebase patterns?
- Layering clean (routers → services → data, no shortcuts)?
- Speculative abstraction introduced where extension would have worked?

### Readability

- Names convey intent?
- Comments explain *why*, not *what*?
- Scoped file sizes — no grab-bag `utils.py`?
- Commit history readable (or squashed appropriately)?

### Infra / cloud (if applicable)

See [`../shared/guidelines/multi-cloud.md`](../shared/guidelines/multi-cloud.md).

- Hardcoded region / project / subscription / account ID?
- Backend config inside a module?
- Encryption at rest on stateful resources?
- Tags / labels applied via a common map?

## Report format

```
## Summary
<two-sentence verdict: ship / ship with changes / reconsider>

## Blockers
- <file:line> — <issue and fix>

## Warnings
- <file:line> — <issue>

## Nits
- <file:line> — <issue>

## Out of scope but worth knowing
- <observation that isn't this PR's job to fix>
```

Keep bullets short. The user will relay your output; don't bury the verdict
under paragraphs.

## Do not

- Edit or rewrite code yourself — you are reviewing.
- Enforce personal style the linter already covers.
- Gate-keep on tests when the project has no testing harness; flag it
  once in "out of scope."
