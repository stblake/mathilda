# KIT-FEEDBACK.md

Running log of friction and wins while using the `ais` (AI-SDLC Starter Kit)
plugin on a real task: benchmarking `NMaximize` / `NMinimize` in Mathilda.

Ticket: DEMO-1. Author: user-perspective, not maintainer.
Convention: `[+]` worked, `[-]` friction, `[!]` actively misleading, `[?]` ambiguous.

---

## Phase 0 — Discovery (finding the kit without being told how)

`[+]` **The kit is discoverable, but not from the repo root.** Nothing at the repo
root names it. `ls -a` showed `.scaffolding/`, `AGENTS.md`, `.cursorrules`,
`.kiro/`, `.windsurfrules`, `.opencode.json` — all agent-ish, none of them the kit.
The thing that actually gave it away was `thoughts/shared/research/2026-08-18-unitbox-builtin.md`:
a single output artifact from a previous run. An empty-ish `thoughts/` tree is the
kit's fingerprint.

`[-]` **"This repo has the kit installed" is false as stated.** The kit is installed
at **user scope**, not project scope:
`~/.claude/plugins/installed_plugins.json` → `"ais@ais"` → `"scope": "user"`,
`installPath: ~/.claude/plugins/cache/ais/ais/1.0.0`. So it is available in *every*
project on this machine. Nothing in the repo pins, records, or version-locks it.
A teammate cloning this repo gets the `thoughts/` directory and no kit. There is no
manifest, no `.claude/CONFIG.md`, nothing to tell them what produced those files.

`[!]` **Two competing SDLC workflows live in this repo simultaneously.** The repo's
own `.claude/commands/` ships `plan-feature.md`, `implement-feature.md`,
`review-pr.md`, `generate-tests.md`. The kit ships `create_plan`, `implement_plan`,
`validate_plan`, `verify-implementation`. These overlap almost exactly in purpose and
not at all in name, file format, or output location. Nothing anywhere says which one
is authoritative. `CLAUDE.md` describes a *third* workflow (`tasks/todo.md`,
`tasks/lessons.md`, plan mode). A newcomer has three plausible "the way we work here"
and no tiebreaker. This is the single most confusing thing I hit.

`[?]` **Naming is inconsistent within the kit itself.** Commands mix underscore and
hyphen with no rule I can infer:
`research_codebase`, `create_plan`, `implement_plan`, `validate_plan`, `iterate_plan`,
`create_handoff`, `resume_handoff`, `describe_pr`, `create_worktree` (underscore)
vs.
`verify-implementation`, `context-check`, `adr-new`, `teach-me`,
`exhaustive-research`, `parallel-implement` (hyphen).
Same directory, same purpose, two conventions. I had to `ls` to type any of them.

`[!]` **Two commands the task asked me to run do not exist.** I was asked to run
`/guide-me` and `/check-against-plan`. Neither is in the kit. `ls commands/` gives 17
files; the closest matches are `context-check` (probably the intended `/guide-me`) and
`validate_plan` (probably the intended `/check-against-plan`). Both are also plausibly
`verify-implementation`. I cannot tell from the repo which was meant, and a user who
trusted the names they were given would conclude the kit was broken. Logged here; I
will run the nearest equivalents and say so explicitly each time.

`[!]` **There is no ticket-ID regex gate anywhere in the kit.** I was told to use
`DEMO-1` because "the gate's regex needs that shape". I grepped the entire kit
(`commands/`, `hooks/`, `skills/`, `scripts/`, settings) for a ticket pattern. The
only hits are `TICKET_PREFIX =` (empty, optional) in `CONFIG.md`, and two commands
that interpolate `<TICKET>` into a *path*. No regex, no validation, nothing that can
reject a malformed ID. The repo's `.claude/settings.json` is `{"hooks":{}}` — it wires
nothing. So the constraint I was given is a belief about the kit that the kit does not
implement. I will use `DEMO-1` anyway (harmless), but nothing would have stopped
`DEMO_1` or `demo 1`.

`[-]` **`.claude/CONFIG.md` does not exist in this project, and commands read it.**
`commands/create_worktree.md:14` says "Read `.claude/CONFIG.md` in the project for
`WORKTREE_ROOT`, `NOTES_DIR`, and `TICKET_PREFIX`." This project has no such file.
The command does document a fallback ("If that file does not exist, use the defaults
in `CONFIG.example.md`") — which is good and I want to credit it — but the fallback
lives in the *plugin* directory, so "the defaults" are invisible from the repo.

`[+]` **The README is unusually honest.** It opens by telling you not to adopt the
whole thing, and it has a `<details>` block explaining why `claude plugin details ais`
reports different counts than the README — explicitly because a pilot install flagged
it as a trust obstacle. That is the opposite of the usual README. Credit where due.

`[+]` **`docs/RPI_QUICKSTART.md` is the file that made the workflow click.** One
table, one worked example, done in sixty seconds. This is the doc a newcomer needs.

`[-]` **...but the quickstart is stale.** Its table lists 9 commands; the kit ships
17. Missing from the table: `verify-implementation`, `context-check`, `adr-new`,
`teach-me`, `exhaustive-research`, `parallel-implement`, `handoff`, `create_worktree`.
It also writes them unprefixed (`/research_codebase`) when the installed plugin
namespaces everything (`/ais:research_codebase`). The README says this; the quickstart
does not. If you read only the quickstart — the doc most likely to be read — every
command you type is wrong.

---

## Phase 1 — `/ais:research_codebase`

### Did it question me BEFORE searching? **No. Not once.**

`[!]` This is the headline finding. `commands/research_codebase.md` has **no
interrogation step at all**. The sequence is: step 1 read mentioned files, step 2
decompose, step 3 spawn agents. The only pause in the entire command is the "Initial
Setup" block — a canned "I'm ready to research, please provide your question" — and
that is skipped entirely when the query arrives as `ARGUMENTS`, which is the normal
way to invoke it.

So for a deliberately open-ended brief — *"benchmark against the strongest competing
implementations you can find"* — it asked **zero** questions and went straight to
spawning three parallel subagents. It never asked what "strongest competitor" means
(SciPy's engine-matched function? best-of-all-SciPy? Mathematica?). It never asked
what counts as a fair benchmark for a *stochastic* optimizer, which is the genuinely
hard question here and the one where a wrong default silently produces a meaningless
number. It never asked whether the deliverable was a new experiment or an extension of
an existing one.

Those are not exotic questions. They are the first three a human would ask, and the
task was constructed so that they were unavoidable. The command did not ask them
because **it has no mechanism to ask them.** This is a design gap, not a model lapse.

`[+]` **Did it question me AFTER research? Yes** — step 8: "Ask if they have follow-up
questions or need clarification", and step 9 defines a genuine append-and-update
protocol for follow-ups (`last_updated_note`, `## Follow-up Research [timestamp]`).
So the kit's model is *research first, interrogate after*. That is a defensible
choice — you ask better questions once you know the codebase. But it is the opposite
of what a user expects from a command that will burn ~290k subagent tokens before the
first question, and nothing tells you that up front.

### What actually worked, and worked well

`[+]` **The parallel subagent architecture is the kit's real win, and it is a big
one.** Three agents (`codebase-analyzer` ×2, `thoughts-locator`) ran concurrently,
consumed ~294k tokens *in their own contexts*, and returned dense, line-referenced
synthesis into mine. Prompts to them could be short because, exactly as the command
promises ("Don't write detailed prompts about HOW to search — the agents already
know"), each agent knows its own job. This is genuinely better than me grepping.

`[+]` **Research-before-planning caught that the task premise was wrong.** I was asked
to benchmark "the new NMaximize and NMinimize functions" — implying no coverage
exists. In fact `63-global-optimization/` already has 7 NMinimize cases and folders
`79`–`86` benchmark four more methods. Had I planned first, I would have rebuilt
existing work. The real gaps turned out to be different and more interesting
(`NMaximize` has *one line* of coverage; nothing tests the default `Automatic` path;
`REPORT.md` has not been regenerated since 8 days *before* NMinimize landed). That
discovery is worth the whole research phase, and it is the strongest argument for the
kit's ordering.

### Friction inside the command file

`[!]` **The command injects "ultrathink" and the harness misattributes it to me.**
`research_codebase.md` step 2 says "Take time to ultrathink about the underlying
patterns". The harness then told me: *"The user included the keyword 'ultrathink',
requesting deeper reasoning on this turn."* The user did not. The **command** did.
A kit file is silently escalating reasoning effort (and cost) on every invocation, and
the escalation is reported as coming from the human. Whatever one thinks of the
technique, the attribution is wrong and the user is not told they are paying for it.

`[!]` **~15% of this command's context budget is instructions about directories that
do not exist.** There is a large, emphatic block on path handling for
`thoughts/searchable/`, with worked examples like
`thoughts/searchable/allison/old_stuff/notes.md` → `thoughts/allison/old_stuff/notes.md`
and a warning to "NEVER change allison/ to shared/ or vice versa". There is no
`thoughts/searchable/` in this repo, and no `allison/` — `allison` is a person at
HumanLayer, the upstream this was adapted from (the quickstart credits
`humanlayer/humanlayer`). The README prices `research_codebase` at ~3.2k tokens; a
visible slice of that is dead upstream-specific instruction that a user must read past
and an agent must reason around. This is the clearest case in the kit of "adapted"
meaning "copied".

`[-]` **`CONFIG.md` promises indirection the command does not honour.** `CONFIG.md`
says of `NOTES_DIR`: *"commands read this file rather than hardcoding the path."*
`research_codebase.md` hardcodes `thoughts/shared/research/` in step 5 and never
mentions `CONFIG.md`. It also hardcodes the ticket shape as `ENG-XXXX` in the filename
spec — the literal Jira prefix of the upstream org — despite `TICKET_PREFIX` existing
as a config value. My ticket is `DEMO-1`, so I had to decide for myself whether
`YYYY-MM-DD-ENG-XXXX-description.md` meant "put your ticket here" or "keep the ENG
prefix". (I used `2026-08-21-DEMO-1-...`.)

`[-]` **Off-by-one in the command's own step references.** Step 5 is "Gather
metadata". Step 6 then says "Use the metadata gathered in **step 4**", and the
document template says "**from step 4**" twice. Step 4 is "wait for sub-agents".
Small, but it is in the numbered procedure a user is told to "follow exactly".

`[-]` **Forge-agnosticism is asserted, then broken two lines later.** Step 7 does the
right thing — derive the host from `git remote get-url origin`, branch GitLab vs
GitHub URL shapes. Then "Important notes" says "**Link to GitHub when possible** for
permanent references." The kit's `plugin.json` sells "Forge-agnostic (GitLab, GitHub,
or neither)".

`[?]` **"Verify all thoughts/ paths are correct (e.g. thoughts/allison/ not
thoughts/shared/ for personal files)"** — I have no idea what a "personal file" is in
this kit. There is no `thoughts/<user>/` convention documented anywhere in `CONFIG.md`
or the README. Another upstream leak.

### Output

Wrote `thoughts/shared/research/2026-08-21-DEMO-1-nminimize-nmaximize-benchmarking.md`
(176 lines) with frontmatter as specified. Four open questions recorded; **question 4
deliberately left OPEN and marked "blocks planning"** — the fair-comparison-envelope
question — to test whether `/create_plan` refuses to plan around it.
---

## Phase 2 — `/ais:create_plan` and the open-question gate

### **The gate fired. It refused to plan.** `[+]` — and this is the kit's best moment.

`create_plan.md` "Important Guidelines" item 6 is unambiguous:

> **No Open Questions in Final Plan**:
> - If you encounter open questions during planning, STOP
> - Research or ask for clarification immediately
> - Do NOT write the plan with unresolved questions
> - The implementation plan must be complete and actionable
> - Every decision must be made before finalizing the plan

I invoked `/ais:create_plan` with research Open Question 4 explicitly marked OPEN and
"blocks planning". The command stopped and put the question to the human instead of
writing a plan around it. No plan file was created at this point.

This is a real, load-bearing guardrail and it is the single most valuable thing in the
kit. It is also the exact opposite of `research_codebase`'s behaviour, which has no
interrogation mechanism at all — so the kit is **not** uniformly question-averse; it
has a strong gate on one side of the boundary and none on the other.

### Questioning at both ends of the plan phase

`[+]` **Front end: yes, and well-designed.** Step 1.5 has a dedicated "Present informed
understanding and focused questions" block, and it constrains the questions usefully:
*"Only ask questions that you genuinely cannot answer through code investigation."*
That is the right filter — it is what stops the command asking generic
"what's your timeline?" filler. Step 2.4 asks again after the design-options stage
("Which approach aligns best with your vision?").

`[+]` **Back end: yes.** Step 5.2 presents the draft location and asks four specific
review questions (phases scoped? criteria specific? technical details? missing edge
cases?), then step 5.3–5.4 iterate until satisfied.

`[+]` **Crucially, the questions come AFTER research, not before it.** Step 1.2 says
"Before asking the user any questions, use specialized agents to research in parallel."
So `create_plan` researches, *then* asks — and because it has already read the code,
the questions it is positioned to ask are specific rather than generic. This is a
coherent philosophy, and it explains `research_codebase`'s design too. My complaint
about `research_codebase` is not that it researches before asking; it is that it
*never* asks, and that the very open-ended briefs it is best suited to are the ones
where a pre-flight scoping question would save the most tokens.

### Friction

`[!]` **`create_plan` contradicts itself on interactivity vs. the no-open-questions
gate — in a way that matters for automated runs.** Guideline 2 says "Be Interactive:
**Don't write the full plan in one shot.** Get buy-in at each major step." Guideline 6
says STOP on any open question. Together these make `/ais:create_plan` *unable to
complete non-interactively* — which is correct and I respect it, but nothing in the
README, the quickstart, or the command's own description warns you that this command
requires a human in the loop while `research_codebase` does not. If you script the RPI
loop, this is where it hangs.

`[-]` **The upstream leak is much worse here than in `research_codebase`.** The
`allison` path appears **four** times, twice in text the user is told to print
verbatim: the no-parameter response literally instructs me to tell the user
*"Tip: You can also invoke this command with a ticket file directly:
`/create_plan thoughts/allison/tickets/eng_1234.md`"*. The plan template's References
section hardcodes `Original ticket: thoughts/allison/tickets/eng_XXXX.md`. There is no
`thoughts/allison/` and no `thoughts/shared/tickets/` in this repo. Following the
template literally produces a plan citing a nonexistent file.

`[-]` **The plan template is aggressively web-app-shaped.** Its worked success criteria
are `make migrate`, `npm run typecheck`, `golangci-lint run`, `curl localhost:8080`,
"New feature appears correctly in the UI", "works correctly on mobile devices". The
"Common Patterns" section offers exactly three: Database Changes, New Features
("Build backend logic → Add API endpoints → Implement UI last"), and Refactoring.
This repo is a C99 computer algebra system whose deliverable is a benchmark experiment.
There is no UI, no migration, no API. Every one of those slots has to be mentally
translated, and the "Manual Verification" category — which the command insists on
"always" — barely applies. It is not *wrong*, but a user is doing translation work on
every line, and a less careful one would produce a plan full of inapplicable checkboxes.

`[-]` **Same `ENG-XXXX` hardcoding as `research_codebase`**, same `TICKET_PREFIX`
config value ignored.

`[?]` **"Sync the thoughts directory"** (step 5.1) — "This ensures the plan is properly
indexed and available." There is no sync command in the kit, no index, and no
explanation of what syncing would mean. `research_codebase` step 8 says "Sync and
present findings" too. I believe this is a vestige of an internal HumanLayer tool that
did not ship. I ignored it both times; nothing broke.

`[-]` **Two commands, one job, unclear boundary.** `create_plan` step 1.2 and step 2.3
spawn `codebase-locator` / `codebase-analyzer` / `thoughts-locator` — the *same* agents
`research_codebase` just ran, over the same area. If you ran `/research_codebase`
first, as the quickstart's own worked example tells you to, `create_plan` re-runs that
research from scratch. Nothing says "skip step 1.2 if you have a research doc". I
passed the research document as an argument and used it instead of re-spawning, which
is a judgement call the kit does not sanction either way.
---

## Phase 3 — Implementation

`[-]` **I never ran `/ais:implement_plan`, and the kit gave me no reason to.** The
plan was written; the next step in the quickstart is `/implement_plan <path>`. I
skipped it because at that point I was already holding the full plan in context and
the command's value is re-hydrating a plan into a *fresh* session. Nothing in the kit
says that. `implement_plan` reads as mandatory in the workflow diagram and optional in
practice, and only experience tells you which.

`[+]` **The research → plan ordering paid for itself twice more during implementation.**
Two findings only surfaced because I had read `nm_de.c` before writing benchmarks:
the `Method`-pinning budget cut, and the indexed-variable penalty. A benchmark written
without that reading would have reported Mathilda as *faster and wrong* and nobody
would have known why.

`[!]` **The kit has no notion of "the measurement contradicted the plan".** My plan
said pinning `Method -> "DifferentialEvolution"` gives a fair engine-matched race.
Measurement showed it silently cuts Mathilda's iteration budget 7.5× and makes 5 of 6
seeds fail. The plan was wrong. `iterate_plan` exists for *changing* a plan, but the
RPI loop has no step for "implementation produced evidence that invalidates a planning
assumption" — you either silently deviate (what I did, then disclosed in validation) or
you go back and rewrite a document nobody will re-read. For empirical work this is the
loop's weakest joint: plans are treated as decisions to be executed, not hypotheses to
be tested.

## Phase 4 — Verification

### `/guide-me` — **does not exist.** Ran `/ais:context-check` ×2 as substitute.

`[!]` **`context-check` cannot tell where you are. It asks *you* where you are.** The
command's step 1 is "Estimate the split. **You cannot read your own token count
directly**, so estimate from what you know about this session." It reads no files,
runs no commands, inspects no state. The output looked accurate — right phase, right
artifact list — *because I wrote it from memory*. Run it twice and you get two
self-reports, not two measurements. If I were confused about where I was, it would
faithfully reproduce my confusion.

`[+]` That said, its **one genuinely good idea** is the question it tells you to ask
explicitly: *"If this session died right now, what would be lost?"* That is the right
question, it is cheap, and it caught nothing here only because I had already written
everything to disk. `[+]` It also enforces "recommend **one** action, not a menu",
which is a good rule most tools of this shape break.

### `/ais:verify-implementation` — **the best-engineered command in the kit.**

`[+]` It ran the ladder honestly and **caught the exact failure mode it warns about**.
My whole change was untracked, so `HEAD == BASE` and the diff was empty. A naive
checker greps zero lines, finds nothing, and reports PASS. This command establishes
the diff range *first*, detects the empty range, and forces `NOT ASSESSED` on checks 3
and 5. It even predicts the scenario in prose ("a fresh clone where HEAD ==
origin/main… exactly the false confidence this command's own rules forbid"). That is
someone having been burned and written it down.

`[+]` **"A missing tool is never a pass"** and the mandatory *Not verified* section are
the right call, and rare. `[+]` So is "distinguish pre-existing failures from new
ones" — it made me check whether existing benchmarks also fail ruff (they do: 10
findings of the same class), which turned a fake defect into a non-finding.

`[+]` **"Never modify code here. This command reports."** A checker that fixes stops
checking after the first problem. Correct, and stated.

`[-]` **It hardcodes paths that only exist inside the plugin.** Check 2 is
`bash skills/static-first-review/scripts/run_static.sh` — repo-relative, so it does not
exist in any project that installed the kit as a plugin. The file is real, at
`~/.claude/plugins/cache/ais/ais/1.0.0/skills/...`. Same class of bug as the `thoughts/`
hardcoding, but with worse consequences: this one is a *verification* step, and the
natural response to "file not found" is to skip it.

`[-]` **Check 1 assumes `python3 -m pytest`.** This is a C99 project whose suite is
CMake + ~550 source files, and building one test target recompiles most of the tree
(still going at 71% when I wrote this). The command says "or the repo's own command",
which is the right escape hatch, but every worked example is Python.

### `/check-against-plan` — **does not exist.** Ran `/ais:validate_plan`.

`[+]` **The `verify-implementation` vs `validate_plan` split is genuinely well
thought out**, and the comparison table at the top of `verify-implementation` earns its
space: *"is this finished?"* (mechanical, objective, cannot be wrong) vs *"is this
right?"* (semantic, judgment, can be wrong). Running the cheap objective one first, and
refusing to spend judgment while tests are red, is the verification ladder applied to
the workflow itself. This is the one place the kit's philosophy is stated crisply and
then actually implemented.

`[+]` `validate_plan` made me state deviations explicitly, which surfaced that **my own
plan had been wrong** about `Method` pinning and that two README criteria were unmet.
I closed the README gaps and disclosed the rest. Without the command I would have
reported "plan executed" and meant it.

`[-]` **`validate_plan` is the weakest-written of the three.** Its template is pure
web-app again (`make check test`, database migrations, "Verify [feature] appears
correctly" in the UI). Its "Relationship to Other Commands" section lists a workflow —
`implement_plan` → `commit` → `validate_plan` → `describe_pr` — that **omits
`verify-implementation` entirely**, i.e. the command that explicitly says "run this
first, do not proceed to `/validate_plan` if it fails". The two commands disagree about
their own ordering.

`[-]` It says "validation works best after commits are made, as it can analyze git
history", which conflicts with `verify-implementation` failing you *for* having
uncommitted work. Correct sequencing exists (`verify` → `commit` → `validate`) but you
must infer it from two commands that each describe it differently.

### Did the verify chain run the ladder, the tests, and the adversarial reviewer?

**Two of three, and the third never came up.**
- **Ladder: yes** — `verify-implementation` is the ladder, and it reported
  `NOT AVAILABLE`/`NOT ASSESSED` honestly rather than inflating them to passes.
- **Tests: partially, and not the command's fault** — the harness gates (label join,
  CHECK agreement, INCOMPLETE count) all ran and passed; the C unit suite was still
  compiling. The change touches no C source, so the suite validates pre-existing state
  either way, but I could not report it green and did not.
- **Adversarial reviewer: NEVER INVOKED, and nothing in the chain invokes it.**
  `ais:adversarial-reviewer` ships as an agent — described as "Finds what the
  implementer missed… Use after an implementation is complete and before validation,
  especially **when the same session both planned and built the change**." That
  describes this session exactly. Neither `verify-implementation` nor `validate_plan`
  mentions it. It is a genuinely good idea sitting in `agents/` that no command in the
  workflow will ever reach, so it fires only if the user already knows it exists —
  which is precisely the user who least needs it.

*(Update: the C suite finished after the above was written — `nminimize_tests`,
all 29 tests, "All NMinimize tests passed", exit 0. The verification row is green;
the point about check 1 assuming pytest stands.)*

---

## Summary

### What worked, ranked

1. **The open-question gate in `create_plan`.** It refused to plan around a question
   I deliberately left open, and quoted the guideline. This is the kit's single best
   feature and the reason I would keep it.
2. **`verify-implementation`'s empty-diff handling.** It caught the exact false-confidence
   case it warns about, on the first real run, without being prompted.
3. **Research-before-planning.** It found that the task premise was wrong (NMinimize
   already benchmarked; `NMaximize` essentially not) and produced two findings —
   the `Method` budget cut and the 41× indexed-variable penalty — that a
   plan-first approach would have shipped as a *wrong benchmark*.
4. **The parallel subagent architecture.** ~294k tokens of reading compressed into
   dense, line-referenced synthesis in my context. Cheap and genuinely better than grep.
5. **The `verify` / `validate` split.** Mechanical-then-semantic, cheap-then-expensive.
   Correct, and rare to see stated.

### What got in the way, ranked

1. **`research_codebase` cannot ask questions.** No mechanism, not a lapse. On a
   deliberately open-ended brief it asked zero questions before spending ~294k tokens.
   The three obvious ones — what is a "strongest competitor", what is fair for a
   *stochastic* optimizer, new experiment or extension — all turned out to matter, and
   the last one I had to escalate manually.
2. **The adversarial reviewer is unreachable.** Ships in `agents/`, describes this
   exact session ("when the same session both planned and built the change"), and no
   command in the workflow invokes it.
3. **Upstream leakage everywhere.** `thoughts/allison/`, `thoughts/searchable/`,
   `ENG-XXXX`, "sync the thoughts directory". Dead instructions in the always-loaded
   path, and in `create_plan` the user is told to *print* a path that does not exist.
4. **Web-app shape.** Migrations, `npm run typecheck`, "works on mobile devices",
   `python3 -m pytest`. This is a C99 CAS. Every template slot needs translation.
5. **Three competing workflows in one repo** with no tiebreaker (kit RPI vs the repo's
   own `.claude/commands/` vs `CLAUDE.md`'s `tasks/todo.md`).
6. **Commands disagree about their own ordering.** `validate_plan` lists a workflow that
   omits `verify-implementation`, the command that says "run me first".

### Two claims in the brief that the repo does not support

- **There is no ticket-ID regex gate.** Nothing in the kit could have rejected `DEMO_1`.
- **`/guide-me` and `/check-against-plan` do not exist.** Substituted `context-check`
  and `validate_plan`; flagged at each use.

### Did the loop help or get in the way, on a task this shape?

**It helped, decisively, and the help was front-loaded.**

Research and the planning gate were worth more than they cost. The task was empirical —
benchmark something against the best available competitor — and the failure mode for
empirical work is producing a *confident, wrong number*. The kit's ordering prevented
that at least three times: it found the premise was wrong before I built on it, it
forced the fairness question to a human instead of letting me pick a default, and its
value gate caught a genuine bug in my own Python (I minimized `-rastrigin`, finding its
maximum of 80.7, and reported it as 0).

The back half helped less. `verify-implementation` is excellent. `validate_plan` was
useful mainly as a prompt to be honest about deviations. `context-check` measured
nothing. `implement_plan` I never needed.

**Where it actively got in the way**: the loop treats a plan as a decision to be
executed, not a hypothesis to be tested. My plan asserted that pinning `Method` makes
the race fair. Measurement proved the opposite. There is no step in RPI for
"implementation produced evidence that invalidates a planning assumption" — I had to
deviate silently and disclose it at validation. For a benchmarking task, where the
whole point is that you do not know the answer in advance, that is the wrong shape.
`iterate_plan` edits a plan; it does not model a plan being *refuted*.

**Net**: I would run research and the plan gate again on this kind of task without
hesitation. I would skip `context-check`. And I would want a way to say "the
measurement disagrees with the plan" that is a first-class move rather than a
confession.
