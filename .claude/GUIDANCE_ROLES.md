# Guidance roles

Some skills ask for org-supplied guidance documents by name — a security baseline, an
engineering-principles file, PR-review conventions — because that content is genuinely
yours to write, not something a generic kit can ship (see `skills/CLAUDE.md`'s Purity
rule: no client, product, or team-specific content in shipped files).

Skills never hardcode a path to your document. They read this file for a named **role**
and resolve it here instead. Fill in the ones your team has; leave the rest
`not-configured` — a skill whose role is unconfigured says so out loud and proceeds on
general practice. It never fails silently and never guesses at a path that doesn't exist.

A role's value is one of three things:

- a path in this repo, e.g. `docs/security-baseline.md`
- a path in a submodule, e.g. `vendor/org-standards/security-baseline.md`
- the literal `not-configured`

Two things to know before filling this in, both found by an audit on 2026-08-29:

- **Only three roles are read by anything this kit ships**: `code-review-checklist`
  (`static-first-review`, `adversarial-reviewer`), `pr-template` (`/describe-pr`), and
  `architecture-guidance` (`/create-plan`).
  Every other role below is *offered* by `/onboard` and read by nothing yet. That is not
  a bug you need to work around — the roles exist so your own skills can resolve your own
  documents by name — but configuring one does nothing until a skill declares
  `GUIDANCE ROLE: <name>` in its body. Each entry says which case it is.
- Earlier versions of this file credited six skills that do not ship
  (`terraform-module-workflow`, `react-feature-delivery`, `database-migration`,
  `analyze-enterprise-documents`, `test-generation`, `design-review`). Those attributions
  were wrong and are gone. `design-system` and `brand-guidelines` were removed for the
  same reason: nothing offered or read them.

## security-baseline

```
security-baseline = not-configured
```

Your organization's security conventions for code that handles auth, secrets, user input,
or infrastructure.

**Read by:** nothing in this kit yet — `/onboard` offers it and your own skills can read it, but no shipped skill does. Configuring it changes nothing until something declares `GUIDANCE ROLE: security-baseline`.

## multi-cloud-guidelines

```
multi-cloud-guidelines = not-configured
```

Your organization's conventions for resources spanning more than one cloud provider (naming,
tagging, cross-account access).

**Read by:** nothing in this kit yet — `/onboard` offers it and your own skills can read it, but no shipped skill does. Configuring it changes nothing until something declares `GUIDANCE ROLE: multi-cloud-guidelines`.

## documents-and-data-guidelines

```
documents-and-data-guidelines = not-configured
```

Your organization's rules for handling extracted document content — PII, retention, what may
leave the building.

**Read by:** nothing in this kit yet — `/onboard` offers it and your own skills can read it, but no shipped skill does. Configuring it changes nothing until something declares `GUIDANCE ROLE: documents-and-data-guidelines`.

## pr-and-review-guidelines

```
pr-and-review-guidelines = not-configured
```

Your organization's PR and code-review conventions (what a reviewer checks, what blocks a
merge).

**Read by:** nothing in this kit yet — `/onboard` offers it and your own skills can read it, but no shipped skill does. Configuring it changes nothing until something declares `GUIDANCE ROLE: pr-and-review-guidelines`.

## agents-guide

```
agents-guide = AGENTS.md
```

An `AGENTS.md`-style file at your own repo's root — engineering principles specific to that
codebase (extend-before-create, no speculative abstractions, house patterns). This is a
per-consuming-repo file, not a kit-wide one — most teams that have one keep it at their own
repo root and would set this to `AGENTS.md`.

**Read by:** nothing in this kit yet — `/onboard` offers it and your own skills can read it, but no shipped skill does. Configuring it changes nothing until something declares `GUIDANCE ROLE: agents-guide`.

## pr-template

```
pr-template = unset
```

Not actually a gap: `/describe-pr` already discovers a PR description template on its own
(`${NOTES_DIR}/shared/pr_description.md`, then your forge's own convention, then a root
`PULL_REQUEST_TEMPLATE.md`, then a built-in fallback) and treats a missing template as
normal, not an error. This role exists only so the Setup Checklist below has a name for that
lookup — setting it to anything has no effect on `/describe-pr`'s behavior.

**Read by:** `/describe-pr`.

## code-review-checklist

```
code-review-checklist = not-configured
```

Your team's own list of failure patterns from past incidents or house-specific gotchas — the
things a generic review rubric won't catch because they're particular to this codebase or
domain.

**Read by:** `static-first-review`, `adversarial-reviewer`.

## architecture-guidance

```
architecture-guidance = docs/design
```

Where your organization's architecture guidance lives — a docs directory, an ADR set, a
vendored standards repo, or a separate repository entirely. A path or a repo URL both work.

This is the pointer to the asset `OPEN-QUESTIONS.md` says the kit cannot ship. The gap named
there is "PRD to architecture": doing it properly needs a repo map plus a product ontology,
and that is per-organization by nature. The kit cannot supply your architecture guidance. It
can stop behaving as though it does not exist.

**Read by:** `/create-plan`, in its Architecture Impact section — a non-empty line there is a
routing signal, and with this role set the plan cites what your org already decided instead
of only saying that a reviewer is needed.

## context-knowledge-graph

```
context-knowledge-graph = not-configured
```

The name of an MCP server, or a command, that answers structural questions about your
codebase: who calls this, what breaks if it changes, which tests cover it. `not-configured`
is the honest answer for most teams and everything still works.

`docs/practices/mcp-or-skill-or-hook.md` carries a worked example of exactly this shape, as
the one case where an MCP server earns its permanent context cost — a graph is genuinely
stateful, and rebuilding an index per invocation would make it useless. The kit ships no such
server, deliberately, and says so there.

**Read by:** nothing in this kit yet. `/onboard` discovers and records it so the agent knows
the capability exists and can reach for it; a skill of your own reads it by declaring
`GUIDANCE ROLE: context-knowledge-graph`. Setting it changes nothing until something does.
