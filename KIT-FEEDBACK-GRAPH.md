# KIT-FEEDBACK-GRAPH.md

Dogfooding the `ais` (AI-SDLC Starter Kit) plugin on a real Graph/HyperGraph feature in
Mathilda (C99 / CMake, no Python in the target repo). Written live, not reconstructed.

Legend: `[+]` worked, `[-]` friction, `[!]` misleading (looked right, wasn't), `[?]` ambiguous.

**Session facts**
- Plugin: `ais@ais`, version **8.0.0**
- Marketplace source: `https://github.com/ms-bain/ai-sdlc-starterkit.git` (branch `main`)
- `gitCommitSha`: `6a33626d600c28c60c930386b1e9a93212873592`
- Installed: 2026-08-22T23:38:42.964Z, scope `user`

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
