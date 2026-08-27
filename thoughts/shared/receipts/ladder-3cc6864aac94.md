# Verification receipt

- commit: `3cc6864aac940073f86203e741ad66bf3bd2c6fb`
- runner: `skills/verification-ladder/scripts/ladder.py` (rungs are subprocesses; no model ran below the judge rung)

| rung | outcome | command |
|---|---|---|
| static | passed | `make check-c99` |
| typecheck | passed | `SDKROOT=$(xcrun --show-sdk-path) make -j8` |
| unit | failed | `cd tests/build && ctest --output-on-failure` |

**RECOMMENDATION:** BLOCKED — unit failed; no verification claim can be made until it is green
"Review-ready" is never emitted by this tool — that verdict requires the judge rung, which is a human or an LLM, not a subprocess this receipt can vouch for.

A receipt is what makes one run servable to every reviewer: the next reader re-runs the diff since this SHA, not the whole review again.
