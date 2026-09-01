# Verification receipt

- commit: `854419997e4b0d8b35fb401783a17b13db2495c9` **plus uncommitted changes** — this run verified the WORKING TREE, which that commit does not fully contain. Commit and re-run for a receipt that vouches for a SHA.
- runner: `skills/verification-ladder/scripts/ladder.py` (rungs are subprocesses; no model ran below the judge rung)

| rung | outcome | command |
|---|---|---|
| static | unavailable | `—` |
| typecheck | unavailable | `—` |
| unit | passed | `make check-c99 && make check-packed-aware` |
| integration | unavailable | `—` |
| judge | unavailable | `—` |

**RECOMMENDATION:** unit-verified — components are not tested together, and nothing is reviewed yet
Not run: static, typecheck, integration.
"Review-ready" is never emitted by this tool — that verdict requires the judge rung, which is a human or an LLM, not a subprocess this receipt can vouch for.

**Judge rung:** an adversarial review ran and left `thoughts/shared/tickets/RG-2/adversarial.md`. Cited, not adjudicated — this tool records that the review exists and where to read it. It does not parse the review's verdicts, so a citation here never means the findings were resolved.

A receipt is what makes one run servable to every reviewer: the next reader re-runs the diff since this SHA, not the whole review again.
