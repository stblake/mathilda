---
description: Analyze an existing codebase and generate a project-level CLAUDE.md that tells Claude Code what the project is, how it's structured, and how to work effectively within it.
argument-hint: [path-to-project]
allowed-tools: Read, Grep, Glob, Bash
---

Analyze the project at `$ARGUMENTS` (default: current directory) and generate
a `CLAUDE.md` file that gives Claude Code the context it needs to work effectively.

A good project `CLAUDE.md` is the single highest-leverage thing a team can
do when adopting Claude Code. Do this once; every Claude session benefits.

## Step 1 — Inventory the project

```bash
ls -la
cat README.md 2>/dev/null | head -60
cat package.json 2>/dev/null
cat pyproject.toml 2>/dev/null
cat go.mod 2>/dev/null
cat Cargo.toml 2>/dev/null
git log --oneline -5 2>/dev/null
```

Identify:
- **What does this project do?** One sentence.
- **Tech stack:** languages, frameworks, databases, cloud provider.
- **Entry points:** where does the app start? (main.py, cmd/, app.py, server.ts, etc.)
- **Test setup:** pytest, vitest, jest, go test, rspec — what's the command?
- **Lint/format:** ruff, eslint, prettier, golangci-lint — how to run clean?
- **Build/run:** how does a developer start the app locally?

## Step 2 — Map the directory structure

```bash
find . -maxdepth 3 -type d | grep -v ".git\|node_modules\|__pycache__\|.venv\|dist\|build" | head -40
```

Identify:
- Source layout (src/, app/, cmd/, lib/, packages/)
- Where tests live
- Where config lives
- Any unusual or non-obvious structure choices

## Step 3 — Identify conventions

Read 3–5 source files at random to find:
- Naming conventions (snake_case, camelCase, PascalCase, kebab-case)
- File organization patterns (feature folders, layer folders, etc.)
- Import style (absolute vs relative)
- Error handling patterns (raise exceptions, return Result types, etc.)
- Logging setup
- Authentication / secrets approach

## Step 4 — Find the commands

```bash
cat Makefile 2>/dev/null | grep "^[a-z]" | head -20
cat justfile 2>/dev/null | head -20
cat package.json 2>/dev/null | python3 -c "import sys,json; d=json.load(sys.stdin); [print(f'{k}: {v}') for k,v in d.get('scripts',{}).items()]" 2>/dev/null
```

Find the exact commands to run tests, lint, format, and start the app.

## Step 5 — Generate the CLAUDE.md

Write a `CLAUDE.md` to the project root using this structure:

```markdown
# CLAUDE.md

## Project
<One paragraph: what the project is, who uses it, what problem it solves.>

## Tech stack
- **Language/runtime:** <e.g. Python 3.12, Node 20, Go 1.22>
- **Framework:** <e.g. FastAPI, React + Vite, Gin>
- **Database:** <e.g. PostgreSQL via SQLAlchemy, Redis>
- **Cloud:** <e.g. AWS, GCP — or "none / local only">
- **Key libraries:** <the 4-5 most important ones>

## Directory layout
<Brief description of the top-level directories — one line each.
Only include directories that aren't obvious from their name.>

## How to run locally
```bash
<exact command(s)>
```

## How to run tests
```bash
<exact command>
```

## How to lint / format
```bash
<exact command>
```

## Conventions
- <naming: snake_case functions, PascalCase classes, etc.>
- <imports: absolute from src/; no relative beyond one level>
- <error handling: service layer raises exceptions; routers catch>
- <testing: unit tests mock external deps; integration tests use test DB>
- <logging: structlog; never print(); include request_id context>

## What NOT to do
- <things that violate conventions or have caused bugs before>
- <anything the team gets wrong repeatedly>

## Out of scope for Claude
- <things Claude should always ask a human about before doing>
- <production systems, billing, security config, etc.>
```

Do not add sections for things you didn't observe. A short, accurate CLAUDE.md
beats a long, speculative one.

## Step 6 — Summarize what you found

After writing the file, report:
- What you observed (stack, conventions, key patterns).
- What you left blank because you couldn't determine it — and what to ask the team.
- Any immediate risks you noticed (no tests, credentials in code, etc.).

## Related

- Template: [`../templates/CLAUDE.template.md`](../templates/CLAUDE.template.md)
- Guide: [`../CONTRIBUTING.md`](../CONTRIBUTING.md) — how to add skills after onboarding.
