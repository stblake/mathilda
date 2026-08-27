---
type: qualify
ticket: STATS-1
date: 2026-08-27
flow: qrispy
score: 10
tier: LOW
basis: projected
method: synthesized `git diff --numstat` for the projected shape of the change (new stats builtin module ~140 lines, init registration ~14 lines, new C test binary ~220 lines, test build wiring ~8 lines), piped to `skills/pr-risk-triage/scripts/score_mr.py --stdin`
---

# Qualify — STATS-1 (core probability & statistics gap fill)

## Score

**10/100 — LOW** (projected; no diff exists yet — scored before research per QRISPY Q).

Top drivers reported by score_mr.py:
- blast radius (+5): 2 top-level modules touched (src, tests)
- diff size (+5): 382 reviewable lines projected

## What the tier picks downstream

- Research fan-out: locator-class first, small (1-2 Bash-capable agents — kit README
  documents the rg/posix_spawn failure for search-only agents; Bash-capable dispatch is
  the prescribed interim fallback).
- Adversarial review: still runs before commit (the kit is explicit: LOW does not drop
  the adversarial pass — "rungs are the review" only covers what rungs can see). Routed
  lean.
- Verification: full ladder per .claude/VERIFICATION_LADDER.md (real C toolchain),
  baseline reproduced before planning.

## Not assessed (per the skill's own honesty rule)

- Intent/correctness of the projected file list — research may move the real
  registration point and test wiring, which would change paths but not materially the
  score class.
- Numerical semantics risk (Mathematica-compatible edge cases) is invisible to a
  path-based score.
