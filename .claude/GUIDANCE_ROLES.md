# Guidance roles

Some skills ask for org-supplied guidance documents by name — a security baseline, an
engineering-principles file, PR-review conventions — because that content is genuinely
yours to write, not something a generic kit can ship.

Skills never hardcode a path to your document. They read this file for a named **role**
and resolve it here instead. Fill in the ones your team has; leave the rest
`not-configured`.

Written by `/setup-kit` on 2026-08-22. DETECT (`kit-setup/scripts/detect_guidance_roles.py`)
found one hit — root `AGENTS.md` — and no signal for the other five; confirmed with the
human rather than guessed.

## security-baseline

```
security-baseline = not-configured
```

## multi-cloud-guidelines

```
multi-cloud-guidelines = not-configured
```

## documents-and-data-guidelines

```
documents-and-data-guidelines = not-configured
```

## pr-and-review-guidelines

```
pr-and-review-guidelines = not-configured
```

## agents-guide

```
agents-guide = AGENTS.md
```

Root `AGENTS.md` documents the code-review-graph MCP tool usage convention for this repo.

## pr-template

```
pr-template = not-configured
```

Not actually a gap — `/describe-pr` discovers its own template; this role has no effect on
its behavior.

## code-review-checklist

```
code-review-checklist = not-configured
```
