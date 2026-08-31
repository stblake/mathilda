# Kit configuration

One file, read by the RPI commands. Change a value here instead of editing commands.

## Notes directory

```
NOTES_DIR = thoughts/
```

Where research findings, plans, and handoffs are written. The layout underneath:

```
thoughts/
└── shared/
    ├── tickets/         # everything about one ticket lives together
    │   ├── <TICKET-ID>/
    │   │   ├── context.md        # the ticket's own description, if kept in-repo
    │   │   ├── research.md       # /research-codebase output (long doc)
    │   │   ├── research-summary.md
    │   │   ├── plan.md           # /create-plan output (long doc)
    │   │   ├── plan-summary.md
    │   └── _archive/<TICKET-ID>/  # moved here by the archive-tickets skill, never deleted
    ├── subsystems/      # durable per-subsystem architecture docs — one per subsystem slug,
    │   └── <slug>.md    # not per ticket. Declares that subsystem's own skills, conventions,
    │                    # test kit, and interdependencies. Read by /research-codebase and
    │                    # /create-plan for whatever subsystems a ticket's front matter names;
    │                    # scaffolded and looked up by skills/subsystem-docs/. Unlike a ticket
    │                    # folder, this is never archived alongside a closed ticket — it
    │                    # accumulates knowledge and changes only when the subsystem's own
    │                    # conventions genuinely change.
    ├── services/        # service ownership and dependency map — one entry per service, a
    │   └── <slug>.md    # pointer index (owner, what it calls, what calls it, a KB pointer),
    │                    # never a copy of another repo's detail. Category-prefixed
    │                    # (skill--grill-me) for this repo's own skills/commands/hooks;
    │                    # free-form for a real, human-authored cross-repo service. Scaffolded,
    │                    # Drafts are never written here automatically, only via --new or a
    │                    # human edit.
    ├── research/    # /research-codebase output when there's no ticket to key by
    ├── plans/       # /create-plan output when there's no ticket to key by
    ├── handoffs/    # /handoff output
    └── prs/         # /describe-pr output when no forge CLI is available
```

`NOTES_DIR` is the only variable here — everything under it is a fixed, documented subpath,
not a separate config entry. A ticket ID is the key: when a plan or research doc names one,
it's filed under `tickets/<TICKET-ID>/`; without one, it falls back to the flat `research/` or
`plans/` directories above.

`thoughts/` is the upstream HumanLayer name and the kit's default, kept so this material
lines up with the source it derives from. If your team already has a docs convention,
change the value above — commands read this file rather than hardcoding the path.

## Worktree root

```
WORKTREE_ROOT = ../worktrees/
```

Where `/implement-plan` puts isolated implementation checkouts. Relative to the repo root.
Keep it *outside* the repo directory so worktrees never get picked up by the repo's own
tooling, test discovery, or `git status`.

## Forge CLI

```
FORGE_CLI = gh
```

Which command-line tool talks to your code host. `auto` detects what is installed and
authenticated. Set it explicitly to skip detection:

| Value | Host | CLI |
|---|---|---|
| `glab` | GitLab | `glab` |
| `gh` | GitHub | `gh` |
| `none` | anything else | commands generate output to a file for manual posting |

Nothing in this kit *requires* a forge CLI. With `none`, `/describe-pr` still writes a full
description — you paste it in yourself.

## Kit root

```
KIT_ROOT = /Users/67840/sandbox-r21/kit-delivered
```

Where this kit is installed, needed only when the kit is **not** inside this repo — i.e. a
plugin install. `/ais:kit-setup` fills it in. Leave blank when the kit's own files sit in
the project (dogfooding, or install method 2).

Only the `prepare-commit-msg` git hook reads it, and only because a git hook runs outside
Claude Code and cannot ask the harness where the plugin lives. An environment variable
(`AIS_KIT_ROOT`) is still honoured as a fallback, but config is the default for a reason: an
`export` does not survive relaunching Claude, so the hook went silent for a reason that was
invisible from inside the session.

## Commit provenance

```
COMMIT_PROVENANCE = full
```

Whether `/commit` adds trailers pointing at what verified the change. Off by default —
writing into your git history is your decision, not the kit's.

| Value | Trailer | Reader needs |
|---|---|---|
| `off` | none | — |
| `receipt` | `Ais-Receipt: <path>` | the repo checked out |
| `full` | `Ais-Kit`, `Ais-Ticket`, `Ais-Receipt`, `Ais-Rungs`, `Ais-Model` | nothing but `git log` |

The levels differ by **where the data lives**, not how much of it there is. The receipt file
already records the rungs, the model mix, and the SHA — so `receipt` buys every field
indirectly for one line of history, while `full` inlines them so an exporter can read
provenance from `git log --format='%(trailers)'` on a bare clone with no working tree.

Every field is derived from a file: the version from `.claude-plugin/plugin.json`, the rungs
from the receipt's own table, the models from the harness-written session transcript. A field
whose source is missing is **omitted, never guessed** — a provenance trailer that is
confidently wrong is worse than an absent one. `skills/verification-ladder/scripts/provenance.py`
builds the block; `--level` overrides this setting for one run.

**This setting does nothing on its own.** The trailers are appended by a `prepare-commit-msg`
git hook, which `/ais:kit-setup --install-git-hook` links for you. Turning this on without
the hook installed produces no trailers and no error — `kit-setup` reports the hook's status
on every run for exactly that reason.

`Ais-Model` reports a *mix* (`claude-opus-5=581 claude-haiku-4-5=112`) because subagents run
on their own tiers, and an orchestrator-only count would report a haiku-heavy session as pure
frontier usage — the one number `docs/practices/model-routing.md` claims to care about.

**Nothing here emits `Co-Authored-By`, and there is no setting for it.** Provenance is
telemetry; authorship is credit, and it renders as a real contributor. Your own `CLAUDE.md`
and the harness's `CLAUDE_CODE_ATTRIBUTION_HEADER` already decide that — a key here would
be a third place for the same setting to disagree with itself.

## Ticket prefix

```
TICKET_PREFIX =
```

Optional. If your team prefixes branches with a tracker ID (`PROJ-1234`, `ACME-88`), set it
here and the commands will use it in branch names. Leave blank to skip ticket-based naming
entirely — nothing in this kit requires an issue tracker.
