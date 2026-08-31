# Verification receipt

- commit: `7701df5b36cab11000efad1f3b2dff094ff2316f` **plus uncommitted changes** — this run verified the WORKING TREE, which that commit does not fully contain. Commit and re-run for a receipt that vouches for a SHA.
- runner: `skills/verification-ladder/scripts/ladder.py` (rungs are subprocesses; no model ran below the judge rung)

| rung | outcome | command |
|---|---|---|
| static | failed | `ruff check .` |

**RECOMMENDATION:** BLOCKED — static failed; no verification claim can be made until it is green
"Review-ready" is never emitted by this tool — that verdict requires the judge rung, which is a human or an LLM, not a subprocess this receipt can vouch for.

**Judge rung:** an adversarial review ran and left `thoughts/shared/tickets/RG-1/adversarial.md`. Cited, not adjudicated — this tool records that the review exists and where to read it. It does not parse the review's verdicts, so a citation here never means the findings were resolved.

A receipt is what makes one run servable to every reviewer: the next reader re-runs the diff since this SHA, not the whole review again.
