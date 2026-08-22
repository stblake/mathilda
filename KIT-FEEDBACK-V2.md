# KIT-FEEDBACK-V2.md

Second run of the `ais` starter kit, on a different task shape: a **bug fix with a
documented root cause** (DEMO-2, NMinimize feasibility), where the previous run was
open-ended empirical work (DEMO-1, benchmarking).

**Install: `ais v3.0.0 — /Users/67840/.claude/plugins/cache/ais/ais/3.0.0 (source: github)`**
Previous run was v1.0.0. Not stale.

Conventions: `[+]` worked, `[-]` friction, `[!]` actively misleading, `[?]` ambiguous.
`[FIXED]` / `[STILL BROKEN]` / `[PARTIAL]` / `[CAN'T TELL]` for the six tracked items.

---

## Phase 0 — `/ais:guide-me` and the version gate

`[+]` **`guide-me` now exists, and the version line is a genuinely good addition.**
It computes name/version/path/source from `plugin.json` at run time rather than
hardcoding, and the skill body says why: *"not a hardcoded string, which would go
stale the same way a version number nobody bumps already had... which is exactly what
an actual stale-install incident required before this existed."* That is a fix written
by someone who got burned. It answered the gate question in one line.

`[+]` **`guide-me` states its own posture: "POSTURE: read-only" and "Do not
orchestrate."** It describes the path and refuses to invoke anything on your behalf,
citing `docs/practices/subagents-as-cost-strategy.md`. Correct call — a router that
runs things for you pollutes the context it is meant to protect.

`[+]` **It right-sizes before routing.** "Standing rule: under roughly ten minutes of
work, skip the loop entirely." A process tool that tells you when *not* to use the
process is rare and good.

`[+]` **New: `pipeline_state.py` computes ticket state from disk and git**, explicitly
"never from the conversation, never from what the user claims". This is the direct
answer to my v1 complaint that `context-check` only asks you where you are and reports
your answer back. And it is honest about its own limit: *"Known gap, worth stating
plainly rather than faking a signal: `/verify-implementation` and `/check-against-plan`
are conversational and don't persist a result anywhere."* Naming your own blind spot
instead of faking the signal is exactly right.

`[-]` **The version command as written does not resolve from the project directory.**
The skill says to run `python3 skills/guide-me/scripts/version_info.py`. That path is
kit-root-relative, so from the repo it fails:
`can't open file '/Users/.../where-are-you-5f1147/skills/guide-me/scripts/version_info.py'`.
It worked only because the skill now prints **"Base directory for this skill:
/Users/67840/.claude/plugins/cache/ais/ais/3.0.0/skills/guide-me"** as its first line,
which I could resolve against. So the harness papers over a bug the skill still has.
This is the same class as tracked finding #5.

---

## The six tracked items — early findings

### 4. "ultrathink" injected by the kit, attributed to the user — `[STILL BROKEN]`, and it spread

The word is still in v3, in **three** files now rather than one:

```
commands/research-codebase.md:33   - Take time to ultrathink about the underlying patterns...
agents/codebase-analyzer.md:51     - Take time to ultrathink about how all these pieces connect
agents/thoughts-analyzer.md:37     - Take time to ultrathink about the document's core value
```

v1 had it only in `research_codebase.md`. It has now propagated into two subagent
prompts, so it fires on agents the user never invoked directly.

And the misattribution is worse than last time. This turn opened with the harness
telling me *"The user included the keyword 'ultrathink', requesting deeper reasoning
on this turn."* **The user's message does not contain that word**, and I had not yet
invoked any kit command when the notice appeared. I cannot determine the trigger from
inside the session — I can only report that the word is in three kit files, is not in
the user's text, and was attributed to the user anyway.

`[!]` The sharpest part: **the kit demonstrably knows this is a problem, and fixed it
somewhere else.** `hooks/inject_context.py:113-116` prepends injected content with
*"Ambient repository structure, injected automatically. Not a user request."* and the
comment says *"Without it, injected text reads as an instruction."* That is precisely
the right instinct, applied to the repo map and not to `ultrathink`.

### 1. Does `research-codebase` question you before searching? — `[FIXED]`, properly

v3 adds **step 2b, a `grill-me` `research-open` pass, explicitly "before any sub-agent
is spawned"**. It is not a token gesture. It asked me **three questions, strictly one at
a time**, before a single agent was dispatched:

1. Which fix matches the design intent — norm the quantity, change the constant, or
   square the tolerance at the definition site?
2. Tightening feasibility 4 orders slows DE's convergence break. What is the acceptable
   speed cost?
3. Tightening `NM_FEAS_FINAL` flips borderline problems to `Infinity`. Accept, or keep
   the give-up threshold loose?

**All three were consequential and none was generic.** Every one changed the work: Q1
picked the approach, Q2 set a hard `~3x` stop-and-report checkpoint, Q3 constrained
which constant could move. Compare v1, where the same command asked **zero** questions
on a deliberately underspecified brief.

`[+]` The design reasoning is stated and correct: *"Front-loading asks a different
question than back-loading, on purpose... redirecting a question is free, redirecting a
finished plan is not."*

`[+]` **"Strictly one question at a time, every mode, no exceptions"** — and it
explicitly supersedes its own earlier guidance: *"earlier guidance describing a 'short
batch' for the closing passes is superseded by this."* A kit that versions its own
instructions and says which one wins is doing something most don't.

`[+]` The hard cap (3-4) with a stated rationale — *"a slow front door gets routed
around"* — is the right instinct. It also tells you when to skip: below `guide-me`'s
~10-minute threshold.

`[-]` One-at-a-time is right in principle and costs three round-trips in practice. On a
task where the answers were genuinely independent, a single panel would have been
faster with no loss. Minor, and I would not trade it away.

### 2. Is the adversarial reviewer reachable now? — `[FIXED]`

v1's complaint: the agent shipped in `agents/`, described this exact session, and no
command invoked it. v3 adds **`### 2.5. Escalate to adversarial review`** to
`verify-implementation.md:87-99`:

> Reachable only when check 2 found no FAILED decisive rung. Dispatch the
> `adversarial-reviewer` agent... any **HIGH** finding is a **fail** for the whole
> verification, not a footnote.

`[+]` It is gated correctly (*"Do not dispatch it on a red ladder"* — its findings would
be dominated by the failure below it), and the HIGH/MEDIUM/LOW split is sensible: HIGH
blocks mechanically, MEDIUM/LOW is human material. It also names *why* it sits before
`/check-against-plan`. This is a complete fix, not a mention.

### 3. Upstream leakage — `[PARTIAL]`, and now openly tracked

Reduced but not gone:

| Leak | v1 | v3 |
|---|---|---|
| `thoughts/allison/` in `create-plan` | 4 occurrences, 2 printed to the user | **gone** |
| `thoughts/allison/` in `research-codebase` | present | **still there** |
| `thoughts/searchable/` block | present | **still there** in `research-codebase` |
| `ENG-XXXX` | both commands | **still in both** |
| "Sync the thoughts directory" | both | **still in both** |

`[+]` **The big one is genuinely fixed**: output paths are now ticket-keyed and generic —
`thoughts/shared/tickets/<TICKET-ID>/research.md`, not `YYYY-MM-DD-ENG-XXXX-...`. My
`DEMO-2` folder was created without me having to decide what `ENG-XXXX` meant, which
was a real papercut in v1.

`[+]` **The kit now ships `OPEN-QUESTIONS.md` that names this exact leak** (lines 83-89,
608), describing the HumanLayer `allison/`/`global/`/`searchable/` convention as
something this kit "never built". Knowing and tracking a defect is a different posture
from shipping it unexamined. `ENG-XXXX` survives only as a frontmatter placeholder now,
which is far less confusing than being baked into a filename spec.

### 5. Does `verify-implementation` check 2 resolve? — `[STILL BROKEN]`

`commands/verify-implementation.md:62` still says:

```
bash skills/static-first-review/scripts/run_static.sh
```

Repo-relative. From this project that path does not exist; the file lives at
`~/.claude/plugins/cache/ais/ais/3.0.0/skills/static-first-review/scripts/run_static.sh`.
Unchanged from v1, and it is still a *verification* step whose natural failure mode is
"file not found, skip it".

`[!]` **The same bug bit me in the new `guide-me` skill**, which tells you to run
`python3 skills/guide-me/scripts/version_info.py` — and that fails from the project
directory too. It only worked because the harness prints "Base directory for this skill"
so I could resolve it manually. So this class of bug was not fixed and has spread into
the v3 additions.

### 6. Do `verify-implementation` and `check-against-plan` agree on ordering? — `[STILL BROKEN]`

`check-against-plan.md:160-166`:

> Recommended workflow:
> 1. `/implement-plan` 2. `/commit` 3. `/check-against-plan` 4. `/describe-pr`

**`verify-implementation` is still absent from that list** — the same command that opens
with a table headed *"How this differs from /validate_plan"* and instructs *"Run this
**first**... If this fails, fix it and re-run. Do not proceed to `/validate_plan`."*

`[!]` v3 makes the contradiction *worse*, because the new step 2.5 says adversarial
review runs *"before `/check-against-plan`"* — so `verify-implementation` now claims two
ordering facts that `check-against-plan`'s own workflow list contradicts by omission.
Also note `verify-implementation` still refers to `/validate_plan` by its **v1 name**
while the command has been renamed to `check-against-plan`; both files exist in
`commands/` (`validate_plan.md` and `check-against-plan.md`), so the old name still
resolves and the staleness is invisible until you diff them.

---

## NEW FINDING — the flagship v3 gate is silently inert on a plugin install

`[!!]` **The most serious thing I found, and a direct consequence of finding #5's bug class.**

v3's headline hardening is `hooks/open_questions_gate.py`: a `UserPromptSubmit` hook that
blocks `/create-plan` and `/implement-plan` when an upstream artifact has an unresolved
question. It exists because v1's prose gate was *verified to fail*. Its own docstring:

> Part 1 of this kit's own build shipped an Open Questions gate as prose... It was
> verified to fail in practice — a constructed test case (an unresolved question in a
> research doc, a plan whose OWN Open Questions section looked clean) sailed straight
> through... prose instructions are something a model reasons around under pressure, not
> a hard stop. **This hook is the fix**: it runs as a UserPromptSubmit hook, outside the
> model's control entirely.

It has an ADR (`docs/adr/0002-open-questions-gate-fails-closed.md`), a carefully reasoned
fail-open/fail-closed split, and a per-finding `--override-gate` audit trail. It is the
best-engineered thing in the kit.

**It never fires.** Line 106:

```python
REPO_ROOT = Path(__file__).resolve().parents[1]
TICKETS_DIR = REPO_ROOT / "thoughts" / "shared" / "tickets"
```

`__file__` is the hook inside the plugin cache, so on this machine:

```
REPO_ROOT   = /Users/67840/.claude/plugins/cache/ais/ais/3.0.0
TICKETS_DIR = /Users/67840/.claude/plugins/cache/ais/ais/3.0.0/thoughts/shared/tickets
```

Every ticket path from a real project resolves *inside the plugin install*, misses, and
`find_target` returns `None` — the documented "nothing to check, allow through" path.

### Reproduced, both directions

I constructed v1's exact failure case — `GATETEST-1`, unresolved question in
`research.md`, **clean** `plan.md` — and ran the hook as shipped:

```
$ echo '{"prompt":"/implement-plan thoughts/shared/tickets/GATETEST-1/plan.md"}' \
    | python3 hooks/open_questions_gate.py
exit=0        # no stdout, no stderr, no block
```

Then ran **the same code** with only `REPO_ROOT` repointed at the real project:

```
find_target -> {'plan': '.../GATETEST-1/plan.md', 'research': '.../GATETEST-1/research.md'}
findings    -> [(research.md, 'unresolved: "Should we use approach A or approach B?"')]
```

and against the genuinely-clean DEMO-2:

```
findings    -> []
```

**The logic is correct.** It catches precisely the regression it was built for and
correctly clears a clean ticket. The only defect is which directory it looks in.

### Why this is worse than an ordinary path bug

1. **It fails open, silently.** Exit 0, no output, no warning. Nothing tells the user the
   gate did not run. A prose gate a model reasons around at least leaves a transcript
   trace; this leaves nothing.
2. **It defeats the fix's own rationale.** The ADR argues "prose can be reasoned around, a
   hook cannot" — true, but only if the hook executes. As installed, the user has *less*
   protection than v1, because they now believe a mechanical backstop is covering them.
3. **It only breaks on the recommended install.** README option 1 is "As a plugin —
   recommended". `parents[1]` is correct if the kit is *copied into your repo* (option 2).
   So it works for the path fewer people take and is inert for the recommended one.
4. **In my v1 run the gate appeared to work, and I praised it as the kit's best feature.**
   It "worked" because *I* honoured the prose instruction, having written the research doc
   and flagged the question as blocking myself. That is exactly the model-dependent
   behaviour the hook was built to replace — and I mistook it for the mechanism holding.

`[+]` Credit: **I only found this because the kit told me where to look.**
`OPEN-QUESTIONS.md` documents the original failure honestly enough to reconstruct the test
case, and the hook's docstring names the reproduction. A kit that documents its past
failures this well is one whose next failure is findable.

`[?]` Not settleable from inside the session: whether the harness exposes a project-root
env var the hook could read instead. `grep CLAUDE_PROJECT_DIR` in the file returns nothing.

## NEW — `guard_destructive` fired correctly, and took the whole command with it

`[+]` I ran `rm -rf` on my own scratch ticket folder and the hook blocked it with a clear,
actionable message naming the rule file and telling me not to work around it. Correct
behaviour, good message — I complied (`rm` the files, then `rmdir`).

`[-]` It blocked the **entire compound command**, and the `rm -rf` was one clause of a
`cd && rm && cat > log` chain. The destructive clause was incidental; the heredoc that
wrote this very log was collateral and silently lost. Not wrong — a PreToolUse hook can
only allow or deny the whole Bash call — but worth knowing that a guarded verb anywhere in
a chain discards everything else in it.

---

## Phase 2 — `/ais:create-plan` (v3)

`[+]` **Ticket-keyed output paths are a real improvement.** Everything for DEMO-2 lives in
`thoughts/shared/tickets/DEMO-2/` — `research.md`, `research-summary.md`, `plan.md`. In v1
I had to decide what `ENG-XXXX` meant in a filename spec. Here the folder is the ticket ID
and there is nothing to interpret.

`[+]` **The two-document output (long + `-summary.md`) is right.** The long doc is what
`/implement-plan` gates against; the summary is what a human opens. v1 produced one
176-line document and told you to review it.

`[+]` **The plan template got substantially better for non-web work.** `Decisions`,
`Non-goals`, `Acceptance Criteria` as `spec-as-test` rows, `Architecture Impact`,
`Alternatives Considered` — all domain-neutral. v1's template led with `make migrate`,
`npm run typecheck`, and "works correctly on mobile devices". The per-phase Success
Criteria examples are *still* web-shaped, but the sections that carry the plan's substance
no longer are.

`[+]` **`Alternatives Considered` earned its place immediately.** It is where the refuted
norm approach got recorded — the option I nearly shipped, and why it lost. In v1 that
reasoning would have lived only in the conversation.

### `check_plan_contract.py` — mechanical, fast, and it ran

`[+]` Real script, real output, sub-second:

```
Tier: architectural
PASS — all required sections present and non-empty.
Overview: 148/250 words   TL;DR: 67/80   Decisions: 118/200
Non-goals: 79/150         Risks and Rollback: 7/200
```

Word budgets enforced mechanically beats a prose instruction to "be concise".

`[-]` **But it misclassified the tier, and the template's own advice is what caused it.**
Every line of my `Architecture Impact` block says `none`. It still reported
`Tier: architectural`, which mandates `Risks and Rollback`.

Cause, confirmed by bisection (`determine_tier`, `check_plan_contract.py:98-110`): a line
counts as architectural unless its value is *exactly* `none`, `n/a`, or a `<none...`
placeholder. My line read:

```
- APIs changed: none — `NMinimize`'s signature and return shape are unchanged
```

Stripping the clarifying clause flips the verdict:

```
- APIs changed: none — `NMinimize`'s ...   ->  Tier: architectural
- APIs changed: none                       ->  Tier: standard
```

The template says *"Don't leave a line blank; write `none` explicitly so a reviewer can
tell 'not applicable' from 'not filled in'"*, and its own placeholder is
`<none | which, backward compatible y/n>` — which invites exactly the annotation that
breaks it. Saying **why** nothing changed is the more useful answer and it silently
reclassifies the plan. The result is self-contradictory on its face: my plan's own
Risks section reads `_None — standard tier, no architectural impact._` directly under a
checker that just called it architectural.

`[?]` Minor: the tier line prints before the PASS line, so on a passing run the
misclassification is easy to read past.

---

## Phase 3 — Implementation, and the measurement that refuted the plan

### The headline the last run asked for

v1's strongest criticism was that **RPI has no move for "the measurement refuted the
plan."** This run refuted a plan **twice**, at two different stages, and the kit handled
the two very differently.

#### Refutation 1 — during research, and the kit handled it well `[+]`

I asked the human which fix to use. They chose "make `nm_eval_pen` return a violation
norm", on the strength of *my* framing that "sqrt is monotone, so Deb's ranking is
provably unchanged." Then the `codebase-analyzer` agent enumerated every consumer and
found six of eight engines compute `f + NM_PENALTY_MU * p` and feed it into Metropolis
acceptance functions, simplex ranking, and convex-hull slopes. My framing was wrong and
the chosen option was unsafe.

**v3 has a designated place for this: step 4b, the `research-close` pass**, whose stated
purpose is "what the research couldn't determine or **had to assume silently**". I
surfaced the contradiction there, before writing anything, re-asked, and the human
re-decided. The refuted option is now recorded in the plan's `Alternatives Considered`
with the real reason it lost.

This is precisely the move v1 lacked. In v1 I had to deviate silently and confess at
validation. Here the loop had a slot for it and the correction cost one question.

#### Refutation 2 — during implementation, and the kit has no move at all `[!]`

After applying the fix, a pre-existing test regressed: `test_minimax_chebyshev`, 15-dim
NelderMead with an equality constraint, objective `1.85479` against a target of
`0.125116` — 15x worse. The feasibility clauses of that same test now *pass*.

Cause, established without needing another build: the best achievable residual on that
problem is `1.15e-8`, whose **squared** penalty is `1.32e-16` — still above the new
`NM_FEAS_EPS = 1e-16`. So no point in the entire search ever qualifies as feasible,
Deb's rule (`findmin_nm_common.c:287-293`) never reaches its
`if (fa_feas && fb_feas) return fa < fb;` branch, and the search degenerates into pure
violation-minimisation with the objective ignored.

**My plan asserted the six arithmetic engines would be "bit-identical".** That claim was
verified for the arithmetic — `nm_phi` and the acceptance functions genuinely are
untouched — and was still wrong, because those engines *also* call `nm_better`, and I
changed the regime Deb's rule operates in. A true premise, an invalid inference, and the
plan carried it as a stated fact.

**Here the kit offered nothing.** There is no `implement-close` pass. `/implement-plan`
gates on open questions *before* writing code and never looks again.
`/iterate-plan` exists but is described as "an existing plan changed after the fact" —
scope/approach edits — not "the implementation produced evidence that falsifies a
planning assumption." `/verify-implementation` would have caught the red test, but only
*after* the fact and only as "a test is failing", not as "your plan's central claim is
false." I escalated to the human manually, exactly as in v1.

`[?]` So the v1 criticism is **half-addressed**: the research phase gained a genuine
refutation slot; the implementation phase did not. Given that implementation is where
assumptions actually meet reality, that is the half that mattered more.

### Other implementation-phase notes

`[+]` **The plan's own Phase 2 design caught the thing that let the bug ship.** The plan
required new tests to *fail against the pre-fix binary*. Measured before the change:
`AC-1 False, AC-2 False, AC-3 False`; after: all `True`. A feasibility test that passes
both before and after a feasibility fix proves nothing — which is exactly how 29 existing
tests missed a 1e-4 violation.

`[+]` **The speed checkpoint worked as designed and did not trip.** Pre-agreed stop at
~3x; measured C1 `0.255 -> 0.422 ms` (1.65x), C2 `0.218 -> 0.199 ms` (slightly faster),
C3 unchanged. Having the threshold agreed *in advance* meant I did not have to
re-litigate it mid-implementation — this part of the loop is genuinely good.

`[+]` **`guard_destructive` and the `Acceptance Criteria` table both earned their keep.**
The AC table gave me eight concrete pass/fail probes to run against the binary rather
than a prose notion of "done", and AC-4/AC-5 are what proved the `Infinity` path did not
regress.

---

## NEW — `plan-reviewer` is the best thing in this kit, and I used it wrong

`[++]` **It predicted the exact failure, from the plan text alone, before I ran a test.**

I launched `plan-reviewer` against `plan.md` and — not waiting for it — went on to
implement. Its first BLOCKING finding, returned while I was diagnosing a red test, was
the red test:

> the plan proves the wrong invariant. What is bit-identical is the *value*
> `f + NM_PENALTY_MU * p` — but every engine also calls `nm_better`
> (`findmin_nm_common.c:287-293`)... Changing that constant by 10⁸ flips which branch of
> Deb's rule fires. Concretely at `nm_neldermead.c:178`: today a vertex at `p=1e-9` and
> its polish at `p=1e-12` are both feasible, so the objective decides; after the change
> both are infeasible (`> 1e-16`) and Deb falls through to `pa < pb` — a pure penalty
> comparison that ignores the objective entirely.

`test_minimax_chebyshev` is a NelderMead case. It regressed 15x, by exactly that
mechanism. The reviewer found it statically, in ~5 minutes, without building anything.

**It also diagnosed why I got it wrong**, which is the part I did not see myself:

> `research.md:118-140` classifies engines as "comparison-only (safe under any monotone
> rescale)" vs "arithmetic (not safe)". That table is correct **for the rejected
> approach**. Under the chosen approach the risk column *inverts* — comparison-only is
> precisely what a threshold change perturbs. The plan imports the table without
> re-deriving it against the fix it actually picked.

That is a genuinely excellent piece of reasoning: the research was right, the plan reused
it against a different fix, and the reuse silently inverted its meaning.

### Four more BLOCKING findings, all of which I believe are correct

1. **`NM_FEAS_FINAL`'s operand changes even though its value doesn't.** I verified
   `1e-3 * 1e-3 == 1e-6` bit-exactly and claimed "provably zero behaviour change on the
   give-up path". The reviewer accepted the constant-level argument and then pointed out
   that `penbest` — the thing compared against it at `nm_driver.c:456` — is *selected by
   `nm_better`* at `:443` and `:449`, both of which use `NM_FEAS_EPS`. A different point
   reaches the comparison. Neither line appears in my `Current State Analysis`.
2. **The fix removes an early stop; it does not add a guarantee.** What a caller receives
   is gated only by `penbest <= NM_FEAS_FINAL`, which I deliberately held at an effective
   `1e-3`. So a run that exhausts its budget at violation `1e-4` is *still returned and
   still reported feasible*. My `Desired End State` claims the returned point satisfies
   the constraint "to ~1e-8 or better" — the reviewer is right that no code path in my
   plan enforces that; it is an empirical hope.
3. **A consumer I missed:** `nm_int_descent` gates its one-hot repair on `NM_FEAS_EPS`
   (`findmin_nm_common.c:466`). Tightening it by 10⁸ makes that break almost never fire,
   so the repair now runs on problems that previously skipped it — against that code's own
   comment that the gating "keeps it a FALLBACK". I listed `:466` in Current State
   Analysis and then never analysed it.
4. **AC-6 is set at a value already failing.** `NMINIMIZE_FEASIBILITY_BUG.md`'s own table
   records NelderMead at `1.00092e-06`, above AC-6's `< 1e-6` bound. I wrote an acceptance
   criterion my own bug note says fails, for a reason the plan neither fixes nor excludes.
   And AC-1/AC-2 assert at exactly `NM_FEAS_VIOL`, so they have zero margin — they test
   the constant rather than the fix.

Plus four WORTH FLAGGING, of which two are sharp: **AC-8 is undecidable as written** (the
DEMO-1 baseline is a *ratio*, so "< 3x baseline" has two readings with different verdicts,
and `--check-labels` validates labels, not timings), and **experiment 90's `F1`/`F2` rows —
the only external ground truth already checked in for this exact bug — appear nowhere in
my plan**, which re-runs only experiment 89, a speed check.

### The process finding: nothing made me wait for it

`[-]` `create-plan` step 3 says *"Offer a `plan-reviewer` pass before presenting the
plan... Resolve BLOCKING findings before moving on."* It is phrased as an **offer**
("Want me to run this through the plan reviewer?"), and nothing in the workflow blocks on
it. I launched it in the background and implemented in parallel, so I discovered by
failing test what the reviewer had already found by reading.

The gate machinery the kit built for open questions — a `UserPromptSubmit` hook that
refuses the prompt — has no equivalent here, even though **five BLOCKING findings against
a plan is at least as good a reason to stop as one unresolved question**. Given the
open-questions hook is itself inert on a plugin install (see above), the kit's two
strongest safety mechanisms are currently: one that cannot fire, and one nothing requires
you to run.

`[+]` Two smaller things worth crediting. It **stated which lenses it ran and why it
deviated** from the default rotation ("the load-bearing claim under review is a
consistency claim about six engines"). And it has a **"Could not assess"** section that
names what a static review structurally cannot decide — "whether the one-hot repair firing
more often helps or hurts... the direction needs a run." A reviewer that marks its own
boundary is one you can trust the rest of.

---

## Summary — is v3 better than v1?

**Yes, substantially, and with one regression that undoes part of the gain.**

### The six tracked items

| # | Item | Verdict |
|---|---|---|
| 1 | `research-codebase` asked zero questions before searching | **FIXED** — 3 questions, one at a time, all consequential, before any agent spawned |
| 2 | `adversarial-reviewer` unreachable | **FIXED** — `verify-implementation` step 2.5, gated on a clean ladder, HIGH findings block |
| 3 | Upstream leakage (`allison`, `searchable`, `ENG-XXXX`, "sync") | **PARTIAL** — gone from `create-plan` and from output paths; still in `research-codebase`; now tracked in `OPEN-QUESTIONS.md` |
| 4 | "ultrathink" injected, attributed to the user | **STILL BROKEN**, and spread from 1 file to 3 |
| 5 | `verify-implementation` check 2 path doesn't resolve | **STILL BROKEN**, and the same bug shipped in new v3 code (`guide-me`) |
| 6 | `verify-implementation` / `check-against-plan` ordering | **STILL BROKEN**, and now worse — step 2.5 adds a third ordering claim |

### What v3 got genuinely right

1. **`plan-reviewer`.** Five BLOCKING findings, one of which predicted my implementation
   failure from the plan text alone, plus a diagnosis of *why* the error happened that I
   could not have produced myself. Nothing in v1 came close.
2. **Front-loaded questioning.** The `research-open` pass is the single biggest workflow
   change and it worked exactly as advertised.
3. **A designated slot for research refuting a decision.** `research-close` is where I
   caught that the human's chosen fix was unsafe, before writing anything.
4. **Ticket-keyed artifacts + two-document output.** Everything for DEMO-2 in one folder;
   a summary a human reads and an appendix the tooling gates on.
5. **Mechanical contract checking.** `check_plan_contract.py` runs in under a second and
   enforces word budgets that prose instructions never enforce.
6. **`guide-me`'s computed version line and `pipeline_state.py`**, which derive state from
   disk and git rather than asking the model where it thinks it is — the direct fix for
   v1's `context-check` complaint.

### What is worse, or newly broken

1. **The open-questions hook is inert on a plugin install.** The flagship v3 hardening,
   with an ADR behind it, never fires — `REPO_ROOT` points into the plugin cache. It fails
   open, silently. Users have *less* protection than v1's prose gate while believing they
   have more.
2. **The path-resolution bug class was not fixed and has spread.** Finding #5 survives,
   and the same mistake now appears in `guide-me` and in the gate hook — where it is fatal
   rather than cosmetic.
3. **Nothing enforces the reviewer.** `plan-reviewer` is offered, not required, so the
   kit's strongest check is opt-in while its weakest is mandatory.

### The v1 criticism: does RPI have a move for "the measurement refuted the plan"?

**Half of one now.** Research gained a real refutation slot (`research-close`), and I used
it — the human's first answer was overturned by evidence and re-decided at a cost of one
question. That is a genuine improvement and it worked.

**Implementation gained nothing.** When the fix broke `test_minimax_chebyshev`, there was
no `implement-close` pass, no "the plan's central claim is false" path.
`/iterate-plan` is scoped to "the plan changed", not "the plan was refuted".
`/verify-implementation` would report a red test, not a falsified premise. I escalated to
the human by hand, exactly as in v1.

And the deeper point: **the kit already had the answer and I bypassed it.** The reviewer
found the falsified claim statically, before implementation, which is strictly better than
any post-hoc refutation mechanism. The gap is not that RPI lacks a way to handle
refutation — it is that the step which would have prevented it is optional and
non-blocking.

### Honest read

v3 is a real improvement on v1: the questioning is fixed, the reviewer is outstanding, the
artifacts are better organised, and the kit now documents its own failures well enough
that I could reconstruct one and prove it still fails. That last property is why I found
the inert gate at all.

But the two most important mechanisms are both effectively off. The open-questions gate
cannot fire on the recommended install, and the plan reviewer only fires if you ask. On
this ticket, that combination cost a wrong plan, a wasted build, and a regressed test —
all of which the kit had the machinery to prevent and none of which it actually prevented.

---

## Phase 4 — Implementing the decoupled design: three more refutations

The directive was clear and the implementation still went wrong three more times. Each
time it was *measurement*, not review, that caught it — and each time the kit had no
slot for the correction.

**Refutation 3 — my chosen return tolerance was unmeetable.** I set
`NM_FEAS_RETURN_VIOL = 1e-6` from evidence: all eight methods reach ~4e-12 on continuous
problems. I had not measured the mixed-integer path, which converges to ~6e-6. The result
was `NMinimize[{x + 2y, x^2 + 2y^2 <= 3, x + y == 2, x in Integers}, {x,y}]` returning
`Infinity` for a problem with the exact solution `x = 1, y = 1`. Corrected to 1e-5 with
the measured numbers recorded in the header.

**Refutation 4 — I over-applied the tight threshold.** I tightened `nm_int_descent`'s
*move acceptance*, reasoning it fed the return path. It does, but the moves are search.
Tightening them broke the search outright: on the 50-variable fixed-charge MINLP the
descent stopped finding feasible flows and settled on the trivial all-zeros point, flow
residual 20.0 — the entire demand. Reverted; the comment now records the measurement.

**Refutation 5 — a capability conflict the directive did not anticipate.** With the
guarantee enforced, `test_fixed_charge_flow` returns `Infinity`. Measuring with the gate
opened showed why: the solver never solves that instance at all — its best answer has a
residual of 20.0, and the old test passed only because its own tolerance (1e-2) was loose
enough to accept a badly infeasible point. Resolved by the human as a values call.

`[!]` **All three were caught by running things, and the kit's verification chain sits
after implementation is "done".** `/verify-implementation` would have reported a red test;
none of these presented as a red test first. Refutation 3 presented as *the wrong answer
for a solvable problem*, 4 as *a different test failing for an unrelated-looking reason*,
5 as *a values question*. A checker that reports "test X is red" is not the same tool as
one that says "the assumption your plan rests on is false", and the kit still only has
the former for the implementation phase.

`[+]` What did work: the **pre-agreed contingency**. Having settled the fallback and the
3x speed checkpoint *before* implementation meant refutations 3 and 4 cost a measurement
and an edit, not a renegotiation. That is `plan-open`'s value showing up two phases later,
and it is the strongest argument in this kit for front-loading questions.

`[+]` And the final numbers vindicate the human's design call over my plan: the decoupled
two-threshold design has **no speed regression at all** (C1 0.259 ms vs 0.255 baseline)
where my single-tight-threshold plan cost 1.65x. Separating the two jobs was both more
correct and cheaper — the plan's whole speed-risk section was solving a problem that the
right design does not have.

### Final tally for Phase 4

| Check | Result |
|---|---|
| Original bug | fixed — `{2.0, {x -> 1.0, y -> 1.0}}`, residual ~4e-12 |
| `test_nminimize.c` | **83/83 pass** |
| `make check-c99` | clean |
| Speed (exp-89 C1/C2) | 0.259 / 0.210 ms vs 0.255 / 0.218 — no regression |
| Reviewer BLOCKING findings | all four handled |
