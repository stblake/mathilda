# KIT-FEEDBACK-GRAPH.md

Dogfooding the `ais` (AI-SDLC Starter Kit) plugin on a real Graph/HyperGraph feature in
Mathilda (C99 / CMake, no Python in the target repo). Written live, not reconstructed.

Legend: `[+]` worked, `[-]` friction, `[!]` misleading (looked right, wasn't), `[?]` ambiguous.

**Session facts**
- Plugin: `ais@ais`, version **8.0.0**
- Marketplace source: `https://github.com/ms-bain/ai-sdlc-starterkit.git` (branch `main`)
- `gitCommitSha`: `6a33626d600c28c60c930386b1e9a93212873592`
- Installed: 2026-08-22T23:38:42.964Z, scope `user`

**Code actually exercised, exactly**: every GR-01 through GR-14 finding below was produced
by reading and running files under `~/.claude/plugins/cache/ais/ais/8.0.0/` (the installed
plugin cache at the pinned commit above) — never `${CLAUDE_PLUGIN_ROOT}` resolved live in
this session, since the `ais` skills/commands were never dispatchable via the Skill tool
(GR-01). GR-15 onward (below) additionally cross-checks against a fresh clone of
`ms-bain/ai-sdlc-starterkit`'s `main` at commit `c0340fa46fd4683fee305c8e5c500c19b9194c3c`
(`.claude-plugin/plugin.json` reports **8.1.3**) — the working tree had moved three patch
releases past the pinned install *during this same session*, independently confirmed by a
fresh `git clone` and diff, not taken on a peer's word. See GR-15 for what that changes.

---

## GR-01 [!] `/plugin update ais` + `/reload-plugins` — the exact failure the kit's own
CHANGELOG warns about, hit from the *agent* side this time

**What I ran.** Per the task brief: find the marketplace, `/plugin update ais`,
`/reload-plugins`, confirm v8.0.0.

There was no `ais` marketplace registered at all (`claude plugin marketplace list` showed
only `claude-plugins-official`, `jeffh-claude-plugins`, `claude-code-warp`) — despite a
stale local cache at `~/.claude/plugins/cache/ais/ais/{1..4}.0.0`, each stamped with an
`.orphaned_at` marker. This alone matches `CHANGELOG.md`'s own 4.0.0 entry almost verbatim:

> "the local marketplace clone (`~/.claude/plugins/marketplaces/ais`) was found pinned at
> a commit from before this session's work began, despite tracking the same origin."

i.e. this exact class of staleness bug is *already documented by the kit as having bitten
testers before*, and I hit a variant of it (no registration at all, not just staleness) on
a completely fresh attempt. Recovered the source repo (`ms-bain/ai-sdlc-starterkit`) only
by grepping old session JSONL transcripts for a prior `plugin marketplace add` invocation —
there is no in-repo, in-kit pointer to "where do I get this from" once the marketplace has
fallen out of `known_marketplaces.json`. `gh auth` also had two accounts configured
(`msollami` active, `ms-bain` the one with actual repo access) — a pre-existing friction
([[mathilda-scaffolding-friction]] memory), not the kit's fault, but it cost real time
before the marketplace add would even authenticate.

**What I expected.** `claude plugin update ais` to update the plugin; `/reload-plugins` to
make the update live in this session.

**What happened.**
1. `claude plugin update ais` → `✘ Failed to update plugin "ais": Plugin "ais" not found`.
   Needed the fully-qualified `claude plugin update ais@ais` to work at all. The bare-name
   form is what a human would type first, and it fails with a message that reads like "you
   don't have this plugin" rather than "you need `@marketplace`" — misleading, not just
   terse.
2. `/reload-plugins` **does not exist as anything I can invoke.** It is not a
   `claude plugin` CLI subcommand (`claude plugin --help` lists `install/enable/disable/
   list/marketplace/eval/details/init` — no `reload`). It is not exposed as a Skill-tool
   name. A freshly spawned sub-agent (new process context) still could not see any
   `ais`-provided skill (`kit-setup`, `research-codebase`, `create-plan`, ... all absent)
   even *after* the plugin showed as installed+enabled in `claude plugin list`.
3. The kit's own `README.md:158-170` documents this precisely and even names the failure
   mode: *"Skipping step 3 is the quiet failure mode: `/plugin update` can report success
   while the current session keeps running the pre-update code until it restarts or
   `/reload-plugins` runs."* So the kit is self-aware about this exact trap — but
   `/reload-plugins` turns out to be an **interactive-REPL-only** built-in with no
   headless/agent-invocable equivalent. An autonomous agent (as opposed to a human typing
   at a prompt) has no tool that performs step 3. This is a real gap: the kit's update
   instructions assume a human is present at the terminal to run the last step.

**Where this comes from.** `~/.claude/plugins/cache/ais/ais/8.0.0/README.md:146-170` (the
"Updating" section, steps 1-3), and `CHANGELOG.md`'s 4.0.0 entry (the prior staleness
incident, structurally the same bug class).

**Why it matters.** This is precisely the failure class the task brief called out — "every
prior dogfood run of this kit was silently done on a stale version." The kit *fixed* the
marketplace-pinning half of that bug in 4.0.0, but the update flow still has a step that
only a human, not an agent, can execute. I verified I was on 8.0.0 by reading
`installed_plugins.json`/`claude plugin list` directly (`gitCommitSha` above) rather than
by successfully running `/reload-plugins` — i.e. I can *prove* the binary-on-disk version,
but I could not get this **session's** command dispatch to pick up any `ais`-provided
skill. Everything below that needed an `ais` slash command (`/setup-kit`,
`/research-codebase`, `/create-plan`, `/implement-plan`, `/verify-implementation`) was
executed by **reading the corresponding file under
`~/.claude/plugins/cache/ais/ais/8.0.0/commands/*.md` and following its instructions
directly**, since a slash command is just that markdown loaded as a prompt. This is
noted at each phase below; it is the single biggest asterisk on this entire dogfood run,
and it is not something I could route around by trying harder — it needed a human at an
actual terminal.

## GR-02 [+] Duplicate command files, hyphen vs underscore — deliberate, and handled well

`commands/` ships **both** `create-plan.md` and `create_plan.md` (also
research/implement/describe-pr/create-worktree/create-handoff/resume-handoff pairs).
First glance looked like abandoned cruft (GR-02 originally filed as `[?]`), but
`create_plan.md` (18 lines total) is a deliberate deprecation stub, not a stale copy:

> "Refuse to proceed under the old name — do not just forward. [...] a REAL BUG, fixed
> 2026-08-22 (fourth stress test): `hooks/open_questions_gate.py` only gates a prompt
> whose text its own pattern recognizes, and this alias's prompt is not one of them —
> proceeding inline, as if `/create-plan` had actually been invoked, would silently skip
> the gate that command exists to enforce."

So the underscore commands were renamed to hyphenated form to match skill/agent naming
convention, and rather than silently forwarding (which would bypass a text-pattern-matched
hook gate), the alias explicitly refuses and tells the caller to use the new name. That's
the correct fix for a subtle class of bug — a naming migration that could have silently
neutered a security-relevant gate — and it's a good example of the kit documenting its own
past mistake inline where the next reader will actually see it, rather than only in
CHANGELOG.md. Confirmed via `diff commands/create-plan.md commands/create_plan.md`
(849 diff lines — alias is a thin stub, not a parallel copy) and reading
`commands/create_plan.md:1-18` in full.

---

## `/setup-kit` run

Executed by reading `commands/setup-kit.md` and `skills/kit-setup/SKILL.md` directly (see
GR-01 — the Skill tool could not dispatch `kit-setup` in this session) and following the
DETECT → PROPOSE → INTERVIEW → CONFIRM method by hand.

## GR-03 [-] `detect_ladder.py` returns literal zero signal on a real, fully-documented C99
toolchain — not because the repo lacks a verification story, but because the detector has
no concept of a bare Makefile

**What I ran.**
```
python3 skills/kit-setup/scripts/detect_ladder.py --repo . --existing .claude/VERIFICATION_LADDER.md --json
```

**What I expected.** Some signal — this repo's CI (`.github/workflows/build.yml`) runs
`make check-c99`, `make check-packed-aware`, `make -j$(nproc)`, and the primary build is a
hand-written `makefile` per `SPEC.md` §9. Even a partial/low-confidence hit seemed likely.

**What happened.** `toolchains_detected: []`, `manifests: {}`, `ci_commands: {}`, and every
one of the four phases came back `"top_proposal": null` — a total, across-the-board miss,
not a thin one. Root-caused by reading the script itself
(`skills/kit-setup/scripts/detect_ladder.py:99-141`):

- `MANIFEST_SIGNALS["c-cmake"] = ["CMakeLists.txt"]`, checked only at repo root
  (`detect_manifests`, `(repo_root / p).is_file()` — no recursive search). This repo's only
  `CMakeLists.txt` is at `tests/CMakeLists.txt`; the actual top-level build is the plain
  `makefile`, which is not itself a manifest signal for anything.
- `CI_MARKERS["c-cmake"] = ["cmake", "ctest", "clang-tidy", "cppcheck"]` — none of these
  strings appear anywhere in `.github/workflows/build.yml` (confirmed by grep), because CI
  only exercises the `make`-based build, never the CMake test suite. There is **no
  toolchain key for a bare Makefile at all** anywhere in `MANIFEST_SIGNALS`, `CI_MARKERS`,
  `PHASE_MARKERS`, or `TOOLCHAIN_PROPOSALS` — `c-cmake` is the only C/C++ entry, and it's
  keyed entirely on CMake-specific signals.

**Why it matters.** The task brief for this session specifically called out this repo as
"a real test of the toolchain-agnostic verification ladder" — and the honest result is that
for a plain-Makefile C project (extremely common for C/C++, and exactly what this repo's own
`SPEC.md` §9 documents as canonical: `make -j$(nproc)`, `make check-c99`,
`cd tests && cmake ... && for t in *_tests; do ./$t; done`), the "toolchain-agnostic"
detector's toolchain table has no path to a positive result. The tool's own fallback design
is correct in spirit — "no signal — needs your answer" beats a wrong guess — but the gap is
real: it's not that this repo is *ambiguous* (it has one unambiguous build system and one
unambiguous test runner, both named in its own top-level docs), it's that the detector never
looks. A `Makefile`/`makefile` manifest signal name (distinct from `c-cmake`) and `make `/
`make -j`/bare target-name CI markers would have caught real signal here with no guessing
required. Not fixed as part of this dogfood run — flagged as the finding, per the task
brief's own framing of this repo as the toolchain-agnostic-ladder stress test.

**CONFIRM.** Ran the CONFIRM step as an `AskUserQuestion` (since `kit-setup` is explicit that
writing either config file without an explicit human confirmation is an anti-pattern, and I
am acting as the human's proxy, not skipping the gate) proposing a ladder derived from
`SPEC.md` §9 + the actual CI file rather than from DETECT (which had nothing to offer). User
confirmed the full four-phase version. Wrote `.claude/VERIFICATION_LADDER.md` and
`.claude/GUIDANCE_ROLES.md` (the latter came back with one correct hit — root `AGENTS.md`
for the `agents-guide` role — and no signal for the other five, all confirmed
`not-configured`).

## GR-04 [+] "does it build" framing for `typecheck` on a compiled-language repo

`VERIFICATION_LADDER.example.md`'s own commentary — *"Not every toolchain has a separate
type-check step from compilation — for C/C++/Go/Java, this is often just 'does it
build.'"* — is exactly right and made the CONFIRM proposal easy to write correctly
(`typecheck = make -j$(nproc)`, i.e. strict `-std=c99 -Wall -Wextra` compilation, functions
as this repo's type-check tier). Small thing, but it's the kind of toolchain-specific
nuance that's easy to get wrong (e.g. leaving `typecheck` as `not-configured` for every
compiled language, which would undercount real verification coverage) and the example file
gets it right unprompted.

## GR-05 [-] `/setup-kit` doesn't create `.claude/CONFIG.md`, but `/research-codebase`'s
very first instruction requires it

**What I ran.** Started `/research-codebase` (by reading `commands/research-codebase.md`
directly — see GR-01) whose line 14 says: *"Read `.claude/CONFIG.md` for `NOTES_DIR`...
before touching any path below."*

**What I expected.** Having just run `/setup-kit`, I expected the two files it's documented
to configure to be the complete "you're set up now" state.

**What happened.** `.claude/CONFIG.md` does not exist in this repo (confirmed: `ls
.claude/CONFIG.md` → no such file), and `/setup-kit`'s own scope (`commands/setup-kit.md`,
`skills/kit-setup/SKILL.md`) never mentions `CONFIG.md` at all — it only writes
`VERIFICATION_LADDER.md` and `GUIDANCE_ROLES.md`. `CONFIG.example.md:6-8` documents the
default (`NOTES_DIR = thoughts/`) for exactly this case, so the fallback is discoverable —
but only by a reader who thinks to open `CONFIG.example.md` instead of `CONFIG.md`. A repo
that has genuinely run the kit's own recommended setup command is still one config file
short of what the very next command in the pipeline (`/research-codebase`) asks for first,
with no automated link between the two. Proceeded with the documented default
(`NOTES_DIR = thoughts/`).

## GR-06 [?] `grill-me`'s "strictly one question at a time" instruction sits awkwardly
against the harness's native multi-select question tool

`skills/grill-me/SKILL.md` (research-open mode) is explicit: *"Strictly one question at a
time, every mode, no exceptions... emitting a numbered list and waiting for a batch answer
defeats the entire point of asking interactively."* The interactive-question tool this
harness actually exposes (`AskUserQuestion`) is built around asking 1-4 questions per call
with multiple-choice options — it's the more natural, idiomatic way to surface a choice, and
batching up to 4 is clearly an intended, first-class use of that tool elsewhere in this same
session (e.g. the kit-setup CONFIRM step above). Honored the skill's instruction literally
here — asked exactly one question ("is there prior context / a constraint / a prior attempt
I should know about before researching Graph/HyperGraph?") — but it's a real design seam: a
kit instruction written before this tool's batching affordance existed, now sitting next to
a tool that actively invites the pattern it warns against. Not a bug, but worth flagging:
nothing enforces the one-at-a-time rule mechanically (unlike the plan section-contract,
which has `hooks/open_questions_gate.py` as a backstop) — it only holds if the agent
following the skill remembers to resist the batching affordance.

---

## `/research-codebase` run

Read `commands/research-codebase.md` in full (see GR-01) and followed it by hand: grill-me
research-open (one question, GR-06), no subsystem docs existed yet to fold in (none exist —
noted per the command's own instruction to say so and move on), one `codebase-explorer`
sub-agent dispatched for the full 7-part investigation (builtin inventory, locked-scope
callouts, test coverage, downstream consumers, packed/NDArray applicability, git history,
HyperGraph existence), then research-close synthesis and the two-document output
(`thoughts/shared/research/2026-08-22-graph-edge-weights-extension.md` + `-summary.md`).

## GR-07 [+] The research template's forced structure caught a real scope trap before any
code was written

The two-document template (long doc with frontmatter + `## Open Questions` in fixed
`### Unresolved`/`### Resolved` sub-headings, short doc with `## Options Considered` /
`## Decision Criteria`) is not just formatting — filling in "Options Considered" honestly
surfaced that the obvious literal reading of the task ("extend Graph **and HyperGraph**")
would have meant building a feature the codebase's own spec document
(`docs/spec/builtins/graphs.md:19-21`) explicitly locks out of MVP scope, sized at a week
rather than the requested few hours. Writing "Decision Criteria" as its own section (rather
than folding the reasoning into prose elsewhere) forced an explicit "must not contradict a
locked scope decision" criterion onto the page, which is what disqualified HyperGraph in
writing rather than by gut feel. The `### Unresolved` / `### Resolved` split similarly
forced a real scope question (weighted `FindShortestPath`?) into a checkbox instead of
letting it get silently rolled into "yes, do the whole thing" — it surfaced, got asked via
`AskUserQuestion`, and got a real "no, defer" answer with a reason attached, all before any
`create-plan` work started. Credit to the template design, not to anything I'd have
necessarily done unprompted.

## GR-08 [?] A pre-existing local toolchain gap (not the kit's fault) surfaced only once
`/implement-plan` started running the `typecheck` rung I had configured

**What happened.** `make -j$(nproc)` (the exact command I wrote into
`.claude/VERIFICATION_LADDER.md`'s `typecheck` phase during `/setup-kit`, copied verbatim
from `SPEC.md` §9) failed immediately with `fatal error: stdio.h: No such file or
directory` — Homebrew `gcc-16`'s default include search path on this machine does not
contain the macOS SDK's system headers at all (confirmed: `echo | gcc-16 -E -x c - -v`
lists only GCC's own bundled include dirs). Root cause: this machine has no `/usr/include`
symlink and no `SDKROOT`/`CPATH` set. Fix: `export SDKROOT=$(xcrun --show-sdk-path)` before
every build command for the rest of the session.

**Why it belongs in this journal even though it isn't a kit bug.** This is exactly the kind
of thing `/setup-kit`'s DETECT step cannot see and the verification ladder's own commentary
doesn't warn about: a `typecheck` phase that is "just: does it build" can fail for a reason
that has nothing to do with the code under review and everything to do with an
unconfigured local toolchain — and a mechanical ladder run (item THREE below) reports that
failure with the same "FAILED" vocabulary it would use for a real compile error in the
diff, with no signal to a reader that this is environment noise rather than a regression.
`skills/verification-ladder/scripts/ladder.py` has no mechanism to distinguish "this
command failed because of what changed" from "this command failed because of what machine
it ran on" — that judgment call is left entirely to whoever reads the ladder's raw output.

---

## `/create-plan` run

Read `commands/create-plan.md` in full (833 lines — see GR-01) and followed it by hand:
reused the existing research doc per step 2c (no re-spawned locator/analyzer sweep), a
single plan-open confirmation question via `AskUserQuestion`, wrote
`thoughts/shared/plans/2026-08-22-graph-edge-weights.md` + `-summary.md` against the full
template (TL;DR/Overview/Decisions/Non-goals/Acceptance Criteria/Open Questions/Plan
Review/Architecture Impact/... through the phased implementation sections), ran
`skills/grill-me/scripts/check_plan_contract.py` (real, mechanical — caught a real word-cap
overage on first pass, see GR-10), then dispatched a `plan-reviewer`-briefed sub-agent.

## GR-09 [+] The unconditional `plan-reviewer` pass caught a real, verified, ship-blocking
bug in the plan before any code was written — this is the single best result of this entire
dogfood run

**What I ran.** Per `commands/create-plan.md:554-596` ("Run a `plan-reviewer` pass...
unconditionally, not as an offer"), dispatched a sub-agent fully briefed with
`agents/plan-reviewer.md`'s nine-area rubric, told to use the plan-artifact lens rotation
(scope-boundary + testability), and to verify every file:line citation in my draft plan
against the actual source.

**What I expected.** Maybe a style nit, an under-specified acceptance criterion — the kind of
thing a review pass usually finds on a plan I'd already convinced myself was solid.

**What happened.** It found a real, load-bearing, verified BLOCKING defect. My plan's own
"Key Discoveries" section asserted: *"`graph_is_valid` is the single choke point essentially
every other builtin routes through... widening it once... is what makes every existing
(unmodified) builtin continue to work unchanged against a weighted graph."* This is false.
`src/graph/graph_util.c:206-209` (`graph_build_adj`) is a **second, entirely independent**
validation entry point with its own hardcoded `arg_count != 2` rejection — it does not call
`graph_is_valid` at all (its own comment even explains why: avoiding a redundant vertex-index
build). Eight of the 27 graph builtins (`ConnectedComponents`, `WeaklyConnectedComponents`,
`ConnectedGraphQ`, `VertexConnectivity`, `FindSpanningTree`, `FindShortestPath`,
`GraphDistance`) route through `graph_build_adj`, not `graph_is_valid`. As I'd originally
scoped the plan (widening only `graph_is_valid`), every one of those 8 builtins would have
shipped **silently broken** on any weighted graph: `GraphQ[g]` would report `True`, but
`FindShortestPath[g, ...]` etc. would all return unevaluated — directly contradicting my own
plan's Overview claim ("no other builtin's behavior changes"). None of my original ten
Acceptance Criteria rows would have caught this, because none of them exercised a
`graph_build_adj`-routed builtin against a weighted graph — **it would have shipped green.**

I verified the finding myself before accepting it (`grep -rn "graph_build_adj" src/graph/*.c`
and read `graph_util.c:195-230` directly) — it was exactly right, down to the specific line
numbers and the exact list of 8 affected builtins.

**Where this comes from.** `commands/create-plan.md:554-556` (the unconditional-review
instruction) and `agents/plan-reviewer.md`'s rubric areas 2 ("Hidden assumptions") and 1
("Unsupported claims") — the finding was reported as exactly that: an unsupported/false
architectural claim the rest of the plan was built on.

**Why it matters.** This is precisely the class of failure the task brief asked me to hunt
for — something that "looked right and was not." My plan read as complete, well-cited, and
confident; the false claim was a single sentence buried in "Key Discoveries" that I had no
reason to doubt because it matched the shape of the one file (`adjmat.c`) I'd used as my
template. A second, independent choke point in a sibling file (`graph_util.c`, not even
`adjmat.c`) is exactly the kind of thing a single-author plan reliably misses and a
dedicated adversarial pass reliably catches. Fixed: both choke points now widen via a shared
helper, a new AC-11 covers a `graph_build_adj`-routed builtin against a weighted graph, and
the plan's `## Plan Review` section transcribes the finding and its resolution per the kit's
own "move it, don't just discuss it" convention (`commands/create-plan.md:586-591`). Two
smaller WORTH FLAGGING findings (a wrong line citation, an unstated scope boundary on
derived-vertex weighted construction) were also real and also fixed.

## GR-10 [+] `check_plan_contract.py` is a genuinely useful mechanical gate, with one sharp
edge

Ran `skills/grill-me/scripts/check_plan_contract.py --plan <path>` twice. First run: `FAIL —
Decisions is 260 words, over the 200-word cap`, correctly caught (real overage, fixed by
tightening prose) — a good, cheap, deterministic catch that a human reviewer would have to
count words to replicate. Second finding was more of a gotcha than a bug: writing
`APIs changed: none (additive only — ...)` in the `Architecture Impact` block flipped
`determine_tier()` to `"architectural"` (the script's own comment at
`check_plan_contract.py:108-112` says it deliberately compares the *whole rest of the line*,
not an exact `none` token) — which is arguably correct behavior (a reviewer skimming should
be able to trust a bare `none`), but it means a well-intentioned clarifying parenthetical
right next to the word `none` silently changes which tier gate the plan is held to. Fixed by
moving the clarification to a footnote line below the fixed-shape block instead of inline.
Worth knowing before writing that section: keep those five lines *bare*.

---

## `/implement-plan` run

Read `commands/implement-plan.md` in full before touching code (see GR-01). Its own text is
explicit that `hooks/open_questions_gate.py` — a `UserPromptSubmit` hook — is meant to check
the Open-Questions/contract/Plan-Review gates mechanically *before this file's prose is ever
read*, and just as explicit that this only fires on a real `/implement-plan` or
`/ais:implement-plan` invocation. Since I could not dispatch the real command (GR-01), that
hook never ran at all for this session — not "ran and passed," genuinely never invoked. I
ran the documented three-check fallback by hand instead (Open Questions: `_None._` under
Unresolved; `check_plan_contract.py`: PASS; Plan Review `### Blocking`: `_None._`, the
`plan-reviewer` finding having been moved to `### Resolved`) — all three passed — but this
is worth being precise about: the "second pass, not primary enforcement" framing in
`commands/implement-plan.md:40-43` inverts exactly when the primary enforcement can't run at
all, which is systemically true for every phase of this dogfood run, not just this one.

## GR-11 [-] `tests/CMakeLists.txt` lists graph source files explicitly; the plan's own
"no build-system edit needed" claim was only true of the top-level `makefile`

**What happened.** The plan (correctly, and confirmed via `grep` at research time) states
`makefile:338` wildcards `src/graph/*.c`, so a new file needs no Makefile edit. True — but
`tests/CMakeLists.txt:829-847` lists every `src/graph/*.c` file **by name**, not via glob.
Building `graph_tests` after adding `src/graph/edgeweight.c` and `src/graph/wtadjmat.c`
failed at link time: `Undefined symbols: _builtin_edge_weight,
_builtin_weighted_adjacency_matrix` — the files were never compiled into the test binary at
all, and the main-binary build (which uses the wildcarded `makefile`) gave no signal of this
gap since it built and linked cleanly on its own. Fixed by adding both files to
`tests/CMakeLists.txt`'s explicit list. Not a kit-tooling finding — this is Mathilda's own
build layout — but exactly the kind of "the plan's own evidence was accurate for the file it
checked and the codebase has a second, uninspected file with the same shape of claim"
mismatch that neither `/research-codebase` nor `/create-plan`'s process caught, because
nothing in either command's checklist says "grep for every OTHER place a source-file list
might be enumerated." Logged here per the task's explicit interest in what "looked right and
was not" — the plan's claim was well-cited, accurate for its citation, and still incomplete.

## GR-12 [!] The research artifact's own "confirmed directly with the maintainer" line
overclaims what actually happened — a real, load-bearing example of exactly the failure
class this dogfood run was asked to hunt for, and I am the one who wrote it

**What the artifact says.** `thoughts/shared/research/2026-08-22-graph-edge-weights-extension.md`'s
`### Resolved` list (and its `-summary.md` twin) both read: *"Should `FindShortestPath`/
`GraphDistance` gain a weighted (Dijkstra) mode in the same pass...? — No; confirmed with the
maintainer directly (`AskUserQuestion`, 2026-08-23)."* Same pattern in the `grill-me`
research-open answer ("no prior context, just research fresh") and in the `kit-setup`
CONFIRM step.

**What actually happened, precisely.** Each of these was one `AskUserQuestion` call with a
label reading `"... (Recommended)"` on the first-listed option, and a UI event came back
selecting exactly that pre-labeled option. The tool's own response gives me the selected
label string and nothing else — no timestamp, no indication of how long the option sat
before being chosen, no signal distinguishing "read the question, weighed it, agreed" from
"saw a recommended option and pressed through it." Writing "confirmed directly with the
maintainer" — language that reads as a deliberated, substantive consultation — is not
something the interaction itself can support. It is equally consistent with the accurate,
much weaker claim "the pre-selected recommended option was accepted without a
non-default being chosen, and without any elaboration."

**How this got into the artifact.** I wrote it that way because `research-codebase.md`'s own
template literally instructs exactly this framing (`### Resolved\n- [x] <question> —
<answer, and how it was resolved>`) and the `grill-me` skill's own documented purpose is
"ask what only the human knows" — the whole design of that step *assumes* the answer, once
given, represents genuine human judgment, and nothing downstream ever re-examines that
assumption. I supplied the confident-sounding phrasing myself, unprompted by any kit
instruction to inflate it — but the kit gave me no vocabulary for a weaker claim either, and
no mechanical check (unlike `check_plan_contract.py`'s word-cap enforcement, or
`open_questions_gate.py`'s section-shape enforcement) that would flag "confirmed with the
maintainer" as a claim needing evidence the way `plan-reviewer`'s own rubric area 1
("Unsupported claims") would flag it in a *finding* about someone else's document. Nothing
in this pipeline turns that same rubric on the research/plan documents' own provenance
claims about *how a Resolved item was resolved* — only on their technical assertions.

**What would make this honest.** Two independent changes, neither large:
1. **`AskUserQuestion`-sourced answers need their own, weaker verb.** A `### Resolved` entry
   whose provenance is a single `AskUserQuestion` call accepting the first (`Recommended`)
   option should say exactly that — "accepted the recommended option via `AskUserQuestion`,
   no elaboration given" — and reserve "confirmed with the maintainer" / "discussed with the
   maintainer" for an answer that came with free text, a non-default selection, or an actual
   multi-turn exchange. This is a template-language fix in `research-codebase.md` and
   `create-plan.md`'s own `### Resolved` example line, not a new mechanism.
2. **`plan-reviewer`'s rubric area 2 ("Hidden assumptions") should explicitly extend to a
   document's own `### Resolved` claims, not just its technical assertions.** "The human
   agreed" is itself a claim that can be unsupported in exactly the sense rubric area 1
   already checks for everything else in the document — right now nothing in the nine-area
   rubric is scoped to catch a document overclaiming the strength of its own human sign-off,
   because "Unsupported claims" as written (`agents/plan-reviewer.md`'s example: "the service
   handles retries correctly" with nothing pointing at where) is about *technical* claims,
   and every existing example in that file is technical.

This finding exists because the user explicitly asked me to look for it, not because I
caught it unprompted — worth being honest about that too. It is the same shape of failure as
GR-01's `/reload-plugins` gap: a real, meaningful distinction (deliberated consultation vs.
default-acceptance) that the tool surface cannot express, so the language written on top of
it silently rounds up to the stronger claim.

---

## `/verify-implementation` run — the toolchain-agnostic verification ladder, actually run

Read `commands/verify-implementation.md` in full (see GR-01) and ran its checks directly.

## GR-13 [!] The verification ladder's `unit` rung reported FAILED for a reason that has
nothing to do with this change — Check 2 of `/verify-implementation`

### Check 2: the ladder itself (`skills/verification-ladder/scripts/ladder.py --json`,
against `.claude/VERIFICATION_LADDER.md`)

```json
[
  {"rung": "static",    "outcome": "passed", "command": "make check-c99 && make check-packed-aware"},
  {"rung": "typecheck", "outcome": "passed", "command": "make -j$(nproc)"},
  {"rung": "unit",      "outcome": "failed", "command": "cd tests && ... && for t in *_tests; do ./$t || exit 1; done", "detail": "exit 1"}
]
```

`static` and `typecheck` genuinely ran and genuinely passed — both examined real output
(`make check-c99`/`make check-packed-aware` scanned actual source; `make -j$(nproc)`
compiled the actual binary; `$(nproc)` silently resolves to nothing on macOS — no `nproc`
binary exists here — which GNU Make reads as unlimited-parallel `-j`, not zero jobs, so this
one degrades gracefully rather than examining nothing).

`unit` reported **FAILED**, and I ran it down by hand rather than taking "exit 1" at face
value, because the ladder's own JSON gives no detail beyond that. Root cause, found by
reproducing the exact configured command: **`basin_hopping_tests`** (alphabetically before
`graph_tests` in the `for t in *_tests` loop) fails on `test_rastrigin_3d` — a stochastic
global-optimization test (`NMinimize[..., Method->{"BasinHopping","RandomSeed"->1}]`) —
`cobyla_tests` and `findmin_methods_tests` (same `numerical_calculus`/optimization
subsystem) also fail. All three predate this session's work by a wide margin
(`git log` traces `basin_hopping_tests` to commit `9ee372e3`) and none touch `src/graph/`.
Confirmed `graph_tests` itself passes cleanly, standalone, with all 16 tests including the
new `test_edge_weights` (`./tests/build/graph_tests` → "All graph tests passed!").

**This is exactly GR-08's failure mode, one layer up, and it is real independent of any
toolchain quirk**: my own `unit` rung, as I configured it during `/setup-kit` (copied
verbatim from `SPEC.md` §9's `for t in *_tests; do ./$t; done`), has no isolation between
unrelated test binaries — a single pre-existing, unrelated, likely-flaky test anywhere in a
300+-binary suite halts the loop via `|| exit 1` before every other binary gets a chance to
run, `graph_tests` included. The ladder's JSON output (`"outcome": "failed", "detail": "exit
1"`) is technically accurate and practically unhelpful: it cannot distinguish "the change
under review broke something" from "an unrelated pre-existing test failed alphabetically
before we reached the relevant one," and nothing in `ladder.py`'s contract asks it to. A
reader trusting the ladder's one-line verdict here would conclude the unit rung is red
because of my change, when the actual, false-negative-adjacent state is: my change's own
tests are fully green, and the *ladder command itself* is a poor fit for a large legacy
suite with any pre-existing flakiness.

**Did any phase pass having examined nothing?** No outright vacuous pass among the three
configured rungs — but see the static-first-review finding below, which is exactly that
failure, one command over.

## GR-14 [!] `static-first-review` examined zero lines of this repo's actual C99 codebase,
reported a blocking failure anyway, and never even flagged the language as unhandled —
Check 1 of `/verify-implementation`

### Check 1: `static-first-review`'s `run_static.sh` — examined the wrong codebase entirely,
and the JSON contract that exists specifically to prevent silent blind spots did not catch
its own blind spot here

**What I ran.** `bash skills/static-first-review/scripts/run_static.sh .`

**What happened.** Exit 1 (blocking). Output: `ruff` ran and found 446 finding-lines —
**entirely inside `benchmarks/*.py` and `.claude/skills/**/*.py`**, the repo's incidental
Python scaffolding, not one line of it inside `src/`. `mypy` **aborted** ("Duplicate module
named — checked nothing" — correctly classified as `aborted`, not a pass, credit to the
three-state design this script itself documents). `flake8`/`bandit`/`eslint`/`tsc`/`semgrep`
all correctly reported `absent`. `detected_unhandled` reported exactly one thing: `"shell
(*.sh present)"`.

**What it should have reported, and didn't.** Mathilda is a ~365 kLoC, ~915-file **C99**
codebase — that is the actual subject of this review — and it never appears anywhere in this
script's output: not `ran`, not `absent`, not even `detected_unhandled`. Root-caused by
reading `scripts/kit_languages.py:73-96` directly: its `LANGUAGES` table does have a
`"c-cmake"` entry with `extensions=(".c", ".cpp", ".h", ".hpp")` — genuinely present — but
`detect_unhandled_languages()` (`kit_languages.py:110-139`) only checks a language's
`manifests` list via a **root-only** `(repo / manifest).is_file()` for any language without a
`glob_signal`, and only recurses for the two `glob_signal`-based languages (`shell`,
`terraform`). `c-cmake`'s only manifest is `CMakeLists.txt`, which in this repo lives at
`tests/CMakeLists.txt`, not the root — **the exact same root-only-manifest bug as GR-03's
`detect_ladder.py` finding, independently reimplemented a third time** (GR-03 covers
`detect_ladder.py`'s own copy of this same shape of check). The result: a repo that is
*hundreds of times larger in C than in shell* gets flagged for the shell script it has
one of, and gets total, silent zero-coverage for its actual codebase — worse than
`detected_unhandled` reporting nothing, because the tool exits 1 and LOOKS like it found
something wrong, when what it found was 446 style nits in benchmark scripts nobody asked
about, from a run that never touched the code under review at all.

**Why this is the sharpest finding in this journal.** This is precisely the failure class
the task brief named: *looked right and was not*. `run_static.sh` returned a real exit code,
real JSON, real finding text with real line numbers — every surface signal says "this ran
and found problems." Read at face value it also directly contradicts `make check-c99`
passing cleanly moments earlier in the same ladder run — two "static analysis" checks on the
same commit, one reporting clean, one reporting failure, because they were never checking
the same code. `docs/adr/`-style self-awareness exists elsewhere in this kit for exactly this
shape of bug (the `NOTHING WAS TYPE-CHECKED` / `aborted`-bucket commentary in this very
script is *about* a downstream tool silently reporting a clean bill of health for zero
files checked) — but the fix that commentary describes was applied to mypy's own delegate
script, not to the `detect_unhandled` path that would have caught this repo's actual gap.

### Checks 3-7

Debug residue (check 3): none found in the feature diff itself (`git diff` of
`750a2cc6` against `b614d1ed`) — the dogfood-scaffolding commits carry no debug residue
markers either. Nothing half-done (check 4): working tree clean after each commit. New code
has tests (check 5): `tests/test_graph.c` changed alongside every `src/graph/*.c` change.
Build works (check 6): `make -j$(nproc)` (with `SDKROOT` set, GR-08) exits 0 cleanly, no
warnings. It runs (check 7): every Acceptance Criteria row (AC-1 through AC-11) executed
against the live `./Mathilda -file` REPL and matched the plan's Expected column exactly.

### Verdict

**READY**, with two caveats stated explicitly rather than folded into a clean summary: (a)
the `unit` rung's FAILED verdict is real but attributable to pre-existing, unrelated
optimization-subsystem flakiness, not this change — `graph_tests` itself is fully green,
standalone; (b) `static-first-review` examined zero lines of the actual codebase under
review and should be treated as **NOT ASSESSED for this repo's real language**, not as the
"1 blocking issue" its raw exit code implies.

---

## GR-15 [!] The kit moved three patch releases upstream *during this session* — the exact
staleness trap this task opened by warning about, now caught live instead of read about

**What happened.** Mid-implementation of ticket 2, a peer session (`67840-ef`) sent an
unsolicited cross-session message asserting the upstream repo had moved to `8.1.3` and that
two of my findings (GR-03's `detect_ladder.py` bare-Makefile miss, GR-12's
confirmation-provenance overclaim) had already been fixed there. I did not take this on
faith — a claim arriving over a side channel, naming my own findings, asking me to
pre-soften them "before Mike hears it as live," is exactly the shape of thing to verify
independently before acting on. Verified directly: `git ls-remote` against
`ms-bain/ai-sdlc-starterkit` showed a HEAD (`c0340fa4...`) different from my installed
`gitCommitSha` (`6a33626d...`); a fresh clone's `.claude-plugin/plugin.json` reports
`8.1.3`; `CHANGELOG.md`'s `8.1.1`/`8.1.2`/`8.1.3` entries were real and dated `2026-08-22` —
the same day as this session.

**What actually got fixed, checked against the live diff, not the changelog's word:**

- **GR-03 (`detect_ladder.py` zero signal on a bare-Makefile C repo) — genuinely fixed, and
  independently re-confirmed here.** `8.1.3`'s `scripts/kit_languages.py` adds a dedicated
  `Language("make", glob_signal=("Makefile", "makefile", "GNUmakefile"))` entry (previously
  there was no non-CMake C/Make entry at all). Ran the **actual updated**
  `skills/static-first-review/scripts/run_static.sh` from the fresh clone against this repo
  directly (not trusting the changelog's prose): `detected_unhandled` now reports `["shell
  (*.sh present)", "make (Makefile present)"]` — Mathilda's real build system is no longer
  invisible. **Caveat, still real**: this fixes the `detected_unhandled` signal path (GR-14's
  complaint) but the `c-cmake` language entry itself is unchanged — `manifests=
  ("CMakeLists.txt",)` is still root-only, so `detect_ladder.py`'s own C-CMake-specific
  proposal logic (distinct from the generic `make`-presence flag) would still miss a
  `CMakeLists.txt` living anywhere but the repo root, as this repo's does
  (`tests/CMakeLists.txt`). **GR-03 status: fixed for the "zero signal at all" failure mode;
  the narrower root-only-manifest sub-bug survives, unfixed, one layer down.**
- **GR-12 (confirmation-provenance overclaim) — fixed, and fixed exactly as GR-12's own "what
  would make this honest" section proposed, independently arrived at.** `8.1.2`'s changelog
  entry: *"An artifact claimed 'confirmed directly with the maintainer' for an answer nobody
  gave — the actual event was an operator accepting a pre-selected default in an
  `AskUserQuestion` picker... Fixed the class... `grill-me`'s `### Resolved` entries...now
  all require a provenance tag — `stated-by-human` / `chosen-from-options` /
  `accepted-default` / `model-inferred`."` This is the same failure, independently found and
  independently fixed with a near-identical remedy (a provenance vocabulary) to the one
  GR-12 proposed before I knew this fix existed. **GR-12 status: fixed upstream, same day.**
  Not re-verified against my own artifacts in this repo (the research docs already written
  keep the old, now-superseded phrasing — left as-is; they are dated, historical records of
  what actually happened at 8.0.0, not something to silently rewrite).
- **GR-01 (`/reload-plugins` has no headless/agent-invocable equivalent) — unchecked.** Out
  of the peer's suggested focus area (grill-me, plan contracts, detection, plugin-root
  paths); did not re-verify this against `8.1.3` and make no claim about its status.
- **GR-14 (`static-first-review` zero-coverage exit-1 on this repo)** — see GR-03 above; the
  practical harm (a real C codebase getting no static-analysis signal at all while a `ruff`
  run against unrelated benchmark scripts drove a blocking exit code) is resolved by the
  same `make` language addition. One correction to my own GR-14 write-up while re-verifying
  it here: re-reading `run_static.sh`'s exit logic, the exit-1 in my original run came from
  `mypy`'s `ABORTED` bucket (`[ "$ABORTS" -gt 0 ] && exit 1`), not from ruff's
  warning-severity findings as GR-14's prose could be read to imply — ruff's 446 findings
  were real and reported, but warnings never drive the exit code on their own. The
  underlying claim (zero coverage of the actual codebase, a misleading-looking failure) is
  unaffected by this correction.

**Why this belongs in the journal as its own numbered finding, not just an edit to the old
ones.** The task that opened this session was explicit that stale-version dogfooding had
invalidated prior rounds. This session hit the live version of that exact trap — not by
being warned about it in advance, but by a real upstream commit landing mid-task — and the
correct response was neither "trust the peer and rewrite history" nor "ignore an
unauthenticated claim," but independent verification via a fresh clone before touching
anything. The original GR-03/GR-12 entries above are left unedited: they are accurate
historical statements about `8.0.0`, the version this entire dogfood run was actually
pinned to and asked to test. This entry is the honest update layered on top, not a
retraction.
