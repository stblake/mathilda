---
name: bug-triage-and-debugging
description: Use when the user reports a bug, shows a stack trace, says "something's broken", or asks for help debugging unexpected behavior. Covers reproduction, isolation, root cause, fix, and regression test — the full loop. Applies to any language or stack.
---

# bug-triage-and-debugging

The default workflow for turning a bug report into a landed fix. Bias
toward reproducing before theorizing.

## When to invoke

- User reports "X is broken" / "Y doesn't work" / "getting error Z".
- User pastes a stack trace, error log, or failing test output.
- User asks "why is this happening?"

## Method

1. **Establish what "broken" means.** What did the user do? What did they
   see? What did they expect? If any of those three are missing, ask.
2. **Reproduce locally.** Get the failure to happen under your control before
   touching any fix. A bug you can't reproduce is a guess.
3. **Narrow the scope.** What's the smallest input / state that triggers it?
   Bisect recent commits if the bug is new.
4. **Find the root cause, not a symptom.** The fix should answer *why* the
   code was wrong, not just how to avoid the error message.
5. **Write the failing test first.** Then make it pass. This proves you
   understand the bug and creates a regression guard.
6. **Fix minimally.** See engineering principles in root [`../../AGENTS.md`](../../AGENTS.md)
   — no tangential cleanup inside a fix PR.
7. **Verify the fix.** Run the full test suite, not just the new test.
8. **Open the PR** with the reproduction in the description so reviewers
   can see what was broken.

## Where to look (a short checklist)

| Symptom                                | Check                                     |
| -------------------------------------- | ----------------------------------------- |
| `null` / `undefined` / `None` error    | Upstream code that returned it; guard     |
| Off-by-one or boundary condition       | Loop edges, empty collections             |
| Race / non-deterministic               | Async ordering, missing await, shared state |
| "Works locally, breaks in prod"        | Env differences, missing env vars, versions |
| Recently introduced                    | `git log -p -- <file>`; `git bisect`      |
| Silent wrong result                    | Test coverage gap; assertion is weak      |
| Flaky test                             | Shared state between tests, time, network |

## Tools by language

- **Python:** `pdb` / `breakpoint()`; `pytest -x --pdb` on failure;
  `pytest -k <name> -vv`.
- **Node/TS:** `--inspect-brk` + chrome devtools; `debugger;` statement.
- **Go:** `delve` (`dlv debug`), `go test -run TestFoo -v`.
- **Rust:** `lldb` / `rust-gdb`; `cargo test -- --nocapture`.
- **Browser/React:** React DevTools, Network tab, `console.trace()` at
  suspected call sites.

## Guardrails

- **No silencing.** Don't catch-and-swallow the exception to make the
  error go away. See the "Errors are information" principle in root [`../../AGENTS.md`](../../AGENTS.md).
- **No `.skip` on a failing test.** If it's wrong, fix it or delete it —
  don't hide it.
- **Don't `--no-verify`** around failing hooks. Fix the hook's complaint.

## Anti-patterns

- Guessing from the stack trace without running the code.
- Applying a fix before the failing test is written.
- A fix that makes the test pass without explaining why the bug happened.
- "Fixed a null check" — that's a symptom fix; find why something was null.
