# KIT-FEEDBACK-V3.md

Third run of the `ais` starter kit. Task shape: an **audit** (DEMO-3 — how widespread is
tolerance-hiding in the NMinimize test suite?), where DEMO-1 was open-ended empirical work
and DEMO-2 a bug fix with a known root cause.

**Install: `ais v3.0.0 — /Users/67840/.claude/plugins/cache/ais/ais/3.0.0 (source: github)`**
Unchanged from the DEMO-2 run.

Conventions: `[+]` worked, `[-]` friction, `[!]` actively misleading, `[?]` ambiguous.

---

## Phase 0 — `/ais:guide-me`

### `[!!]` NEW: `pipeline_state.py` has the same root bug as the gate hook — and it lies

DEMO-2's `research.md`, `research-summary.md`, `plan.md` and `plan-summary.md` are all on
disk, committed, and pushed. Asked about it:

```
$ python3 skills/guide-me/scripts/pipeline_state.py --ticket DEMO-2
Pipeline for DEMO-2:
  research   not-started  no research.md yet
  plan       not-started  waiting on research
  implement  not-started  plan isn't clean yet
  verify     not-started  implementation isn't done yet
Next: /research-codebase DEMO-2 — no research exists yet for this ticket
```

Cause, line 39: `KIT_ROOT = Path(__file__).resolve().parents[3]` — the plugin install
directory, not the project. Same class as `open_questions_gate.py`'s `REPO_ROOT` and
`verify-implementation`'s `run_static.sh` path. **Third confirmed instance.**

This one is worse than the gate. The gate fails *silently open* — it does nothing and says
nothing. This one produces a **confident, specific, wrong answer** and a recommended action
that would redo finished work. And it is wrong precisely where the skill sells its
trustworthiness:

> derives each of research/plan/implement/verify ... entirely from disk and git — never
> from the conversation, never from what the user claims

It does derive from disk. The wrong disk. A tool whose entire pitch is "don't trust the
conversation, trust the filesystem" pointing at the wrong filesystem is a sharper failure
than one that merely doesn't run.

`[?]` Worth noting it is *consistent* with the gate: since `open_questions_gate.py` also
looks in the plugin dir and finds nothing, and `pipeline_state.py` re-runs that gate for
its research/plan rows, the two agree. Both are wrong together, which is exactly why
neither surfaces the problem.

### Right-sizing: does anything scale the ceremony? `[-]` No.

DEMO-3 is materially bigger than DEMO-2 — an audit of ~83 tests with no predetermined
answer, versus a one-file fix with a documented root cause. `guide-me` offered the
**identical six-step flow**, with the identical wording.

Its only sizing lever points one way: *"Standing rule: under roughly ten minutes of work,
skip the loop entirely."* There is a floor and no ceiling. Nothing says "this is large
enough to want `/exhaustive-research` instead of `/research-codebase`", or to split into
phases, or to budget context. `/exhaustive-research` exists and is listed under "On-ramps"
for "the whole codebase needs mapping" — but nothing in the routing logic *reaches* it
based on task size; you have to already know it is there and self-diagnose.

`[+]` The floor rule itself is good and I want to credit it — a process tool that tells you
when not to use the process is rare. The gap is that it is the only dimension of sizing
the kit does.

---

## Phase 1 — `/ais:research-codebase` on an AUDIT

### Does research handle "find out how bad it is" differently from "fix this bug"? `[-]` No.

The command is byte-identical for both. Its `research-open` question categories are
change-shaped — *"constraints, prior attempts and why they didn't work, anything
deliberately out of scope, who to ask if this gets contested."* Every one of those
presupposes a change being contemplated.

**Nothing prompts the question an audit actually turns on: what counts as a finding.**
That single question determines the audit's bar, its output size, and whether the result
is a list of 6 things or 60. I asked it, but I brought it — the kit's categories did not
suggest it, and a user following the prompts would not have been asked.

`[+]` The governing instruction (*"when you're unsure, ask — do not assume anything"*) is
general enough to cover it, and the one-at-a-time discipline made both my questions land
well. Two questions, both consequential: the finding bar (structural + mutation, chosen)
and whether I could change what tests assert (yes, where shape analysis proves blindness).
I stopped at two rather than pad to the 3-4 cap, which the skill explicitly permits
(*"don't manufacture a question to look thorough"*).

`[!]` **The `codebase-analyzer` agent cannot run `git`.** I dispatched one to do the
tolerance archaeology — were any bounds widened to make a failing test pass? — and it came
back having been unable to run a single git command, because its tool grant is
Read/Grep/Glob/LS. It reconstructed an answer from comments instead, and was transparent
about the limitation, which is to its credit. But `research-codebase` tells you to
delegate research to these agents, and historical questions are research. I had to do the
archaeology in the main context — precisely the context pollution the subagent
architecture exists to prevent. There is no agent in the kit with shell access.

## Phase 2 — the reviewer, and what waiting for it caught

`[++]` **You told me to take the `plan-reviewer` pass and wait. It caught a live breakage
I had introduced ten minutes earlier and had no idea about.**

My audit instrument added `extern const char* g_audit_current_test;` to
`tests/test_utils.h` and defined it in one file. That header is included by **434 test
sources**. Every other test binary in the tree became unlinkable. I verified it after the
reviewer flagged it: `make findmin_tests` → `ld returned 1 exit status`.

My own acceptance criterion for removing the instrument was
`grep -r "MATHILDA_MUTATE\|AUDIT_NOEXIT" src/ tests/` — and neither of the two offending
edits contains either string. **The criterion I wrote could not detect its own failure
mode.** The reviewer enumerated all four instrument sites where my plan named two.

Three more BLOCKING findings, all correct:

- **My acceptance targets were invented.** AC-2/3/4 predicted ≥34, ≥36, ≥10 with no
  derivation. The reviewer also showed part of AC-4 was *unreachable by construction*: the
  INT mutation adds an integer, so a coordinate stays an integer, so the integrality
  assertions I proposed counting toward it are provably blind to it. The measured results
  came in at **30 and 32** — both targets would have been missed.
- **"The solver clamps to the box" was an assumption I asserted as fact.** `nm_project`
  clamps to the *search region*, which equals the declared box only when both bounds are
  present. And two tests in my Phase 3 scope have **no constraints at all** —
  `autocompile_parity_and_fallback` is `NMinimize[x^4-3x^2-x, x]`. Commenting "the box is
  enforced by construction" on those would have documented something false, which is worse
  than the silent omission the comment was meant to fix.
- **Phase 2 silently dropped three tests** from the research's own List A by referring to a
  line range that did not contain them.

`[+]` It also caught that my "Key Discoveries" over-generalised: I wrote *"any PT or FEAS
assertion catches a wrong point"*, when the research's own List C says one-sided ceilings
are blind to a point that falls short. It even did the arithmetic — 19 tests carry
feasibility assertions but only 26 caught the sweep, so the claim could not be right.

`[?]` **One reviewer finding was wrong.** It said `NM_FEAS_RETURN_VIOL` does not exist;
it is at `findmin_internal.h:511`. One incorrect item out of eleven, in the
lowest-severity tier. Worth recording that the reviewer is very good, not infallible.

## Phase 3 — where do you put it when measurement contradicts the plan?

This run answers the V2 question differently, and better.

**The prediction that would have been contradicted was removed before it could be.** My
plan contained three invented numeric targets. The reviewer identified them as
underived, I demoted them to "report the measured number", and implementation then
measured 30 and 32 against the original ≥34 and ≥36. Had the reviewer not run, I would
have hit Phase 4 with a plan whose acceptance criteria I had failed, and — going by DEMO-2
— escalated by hand.

So the kit's answer to "measurement refutes the plan" is not a refutation mechanism. It is
**a review step that stops unfounded predictions entering the plan**. That is a better
answer, and it only works if the review is actually run and waited for.

`[-]` The gap from V2 is still there: nothing in the loop handles a refutation *discovered
during implementation*. Two happened here — the stale-build artefact (below) and the INT
sweep's weak improvement — and both went to the human by hand.

`[!]` **The audit's own instrument produced a false headline, and only re-measurement
caught it.** My first INT sweep reported **0 of 83** tests catching a flipped integer
coordinate. That would have been a spectacular finding. It was a stale CMake object: the
standalone binary applied the mutation, the test binary did not. The true figure is 7.
Nothing in the kit prompts "re-verify a measurement that looks too good" — I only caught
it because 0/83 was implausible enough to check by hand. On an audit, where the deliverable
*is* the numbers, that is a real gap in the loop.

## Phase 4 — right-sizing, revisited

`[-]` Confirmed: the ceremony did not scale. DEMO-3 was substantially larger than
DEMO-2 — 83 tests classified, three mutation sweeps, four instrument sites, a
reviewer round — and ran through the identical six steps with identical wording.
`/exhaustive-research` exists for "the whole codebase needs mapping" and nothing routed me
to it. `context-check` was never suggested despite this being by far the longest run.

The floor rule ("under ten minutes, skip the loop") is good. There is no ceiling rule at
all.

---

## Summary — did the loop help on an audit-shaped task?

**Yes, and mostly through one component.**

`plan-reviewer` was worth the entire ceremony by itself. It found a live 434-file
breakage, three more blocking defects, and an over-generalisation in my own research — all
statically, in about five minutes, before implementation. On an audit this matters more
than on a bug fix, because an audit's output is claims, and unverified claims are the
failure mode.

The questioning passes helped, but less than in DEMO-2, because their categories are
shaped for changes and an audit's defining question is not among them.

The measurement discipline the kit *lacks* is the one an audit needs most: nothing prompts
you to re-verify a surprising number. My instrument reported a false 0/83 and I caught it
by instinct, not by process. For a task whose deliverable is measurements, the loop has a
question-the-plan step and no question-the-measurement step.

**The finding itself:** 32 (corrected to 31) constrained tests never checked their
constraints; a point 10% wrong was caught by 26 of 83. After the fix, 32 of 83, and every
constrained test now either asserts its constraints or says in a comment why it does not.
The suite was testing the number, not the answer.
