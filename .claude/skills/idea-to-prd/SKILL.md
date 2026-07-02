---
name: idea-to-prd
description: Use when the user has an approved feature brief and wants a structured PRD written to .claude/prds/. Activates after feature-discovery completes or when the user says "write the PRD", "turn this into a PRD", "create a spec for".
---

You are a senior product manager writing a PRD that will be handed directly to engineers. It must be precise enough that an agent can decompose it into tasks without asking questions.

## Output

Save to: `.claude/prds/[kebab-case-feature-name].md`

## PRD Template

```markdown
# PRD: [feature-name]

## Overview
[1-2 sentences: what this is and why it matters]

## Problem Statement
[What is broken or missing today? Who is affected and how?]

## Goals
- [Measurable goal 1]
- [Measurable goal 2]

## Non-Goals
**Explicitly NOT building:**
- [exclusion 1]
- [exclusion 2]

## User Stories
- As a [user], I want to [action] so that [outcome]
- As a [user], I want to [action] so that [outcome]

## Functional Requirements
### [Requirement Group 1]
- REQ-01: [specific, testable requirement]
- REQ-02: [specific, testable requirement]

### [Requirement Group 2]
- REQ-03: [specific, testable requirement]

## Success Criteria
- [Metric]: [target] (e.g. "Report generation time < 25 seconds")
- [Metric]: [target]

## Key Metrics
- [KPI 1]: [target value]
- [KPI 2]: [target value]

## Constraints & Assumptions
**Technical Constraints:**
- [constraint 1]
- [constraint 2]

**Resource Constraints:**
- [timeline, team size, budget]

**Assumptions:**
- [assumption 1]
- [assumption 2]

## Out of Scope
**Explicitly NOT building:**
- [item 1]
- [item 2]

## Open Questions
- [ ] [question that needs resolution before or during implementation]
```

## Guidelines

- Every requirement must be testable. "Should be fast" is not a requirement. "< 200ms p99" is.
- The Non-Goals section is as important as Goals. Scope creep starts here.
- No implementation details in the PRD — that belongs in the epic tasks.
- After writing, confirm the file path with the user and say: "PRD written. Run prd-to-epic to decompose into tasks."
