---
description: Scaffold a new skill, command, agent, or Cursor rule with correct frontmatter and location.
argument-hint: <kind> <name> [description]
allowed-tools: Bash, Read, Write, Edit, Glob
---

Scaffold a new AgentKit asset. Parse `$ARGUMENTS` as: `<kind> <name> [description]`.

- `<kind>` = `skill` | `command` | `agent` | `rule` | `guideline`
- `<name>` = kebab-case (`my-skill`, `backend-engineer`, `go-service`)
- `[description]` = optional one-liner; ask for it if missing

## Where to create it

First determine context: are we inside an AgentKit checkout (has `skills/`, `commands/`,
`agents/` at the root), or inside a consuming project with `.claude/` installed?

| Context | Path |
| --- | --- |
| Consuming project (`.claude/` present) | `.claude/skills/<name>/SKILL.md` etc. |
| AgentKit library root | `skills/<name>/SKILL.md` etc. |

Full dispatch table (consuming project):

| kind | Creates |
| --- | --- |
| `skill` | `.claude/skills/<name>/SKILL.md` |
| `command` | `.claude/commands/<name>.md` |
| `agent` | `.claude/agents/<name>.md` |
| `rule` | `.cursor/rules/<NN>-<name>.mdc` (ask for `NN` prefix) |
| `guideline` | `shared/guidelines/<name>.md` (AgentKit only) |

## Before writing

1. **Check for collisions.** If the target path exists, stop — don't overwrite.
2. **Validate the description.** Must be trigger-shaped: "Use when the user asks to…".
   Reject summary-shaped descriptions like "Handles X-related tasks."
3. **For rules**, ask for the `globs` array — don't guess the file pattern.

## Use the template

Find the right template from `templates/` (or from where AgentKit is installed):

- Skill → `SKILL.template.md`
- Command → `command.template.md`
- Agent → `agent.template.md`
- Rule → `rule.template.mdc`

Fill every `<placeholder>`. Mark genuinely optional sections with `<!-- TODO -->`.

## After creating

1. **Print what was created** and what still needs filling in (the `<...>` fields).
2. **Remind the user:**
   - Description must be trigger-shaped or discovery won't work.
   - No project-specific names, credentials, or internal URLs in assets meant to be reused.
   - Reference `shared/guidelines/` for security, multi-cloud, or document rules
     rather than duplicating them inline.
   - Keep the body under ~120 lines — it's loaded into context on every use.
