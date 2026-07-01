---
name: generate-pptx-deliverable
description: Use when the user asks to produce a PowerPoint deck — presentations, status updates, data read-outs, or structured deliverables. Takes a JSON/YAML spec and writes a .pptx with consistent slide layouts (title, section, bullets, two-column, table). Ships a runnable `generate.py`.
---

# generate-pptx-deliverable

Turns a structured slide spec into a PowerPoint file. Deck generation becomes
a one-command step so findings flow out of analysis into deliverables without
hand-formatting.

## Use it directly

```bash
pip install -r requirements.txt
python generate.py deck.yaml --output final.pptx
```

Minimal spec (YAML):

```yaml
title: Q3 Market Assessment
subtitle: Preliminary findings
author: Your Team Name
slides:
  - type: title
    title: Q3 Market Assessment
    subtitle: Working session — 2026-04-23
  - type: section
    title: 1. Market size and growth
  - type: bullets
    title: Key findings
    bullets:
      - Addressable market grew ~8% YoY
      - Three incumbents hold ~65% share
      - Adjacent white-space in mid-market (<$50M)
  - type: two_column
    title: Supply vs. demand
    left_title: Supply side
    left_bullets: [Fragmentation, Capacity constraints, M&A activity]
    right_title: Demand side
    right_bullets: [Enterprise shift, Regulatory tailwinds]
  - type: table
    title: Competitor positioning
    headers: [Player, Share, Trajectory]
    rows:
      - [Incumbent A, "30%", Stable]
      - [Incumbent B, "22%", Gaining]
      - [Incumbent C, "13%", Declining]
```

Or equivalently JSON. Spec schema is documented at the top of
[`generate.py`](generate.py).

## When to invoke

- User asks to "make a deck", "draft slides for X", or "produce a PPT from this analysis".
- Another skill produces structured findings (e.g.
  [`analyze-enterprise-documents`](../analyze-enterprise-documents/SKILL.md))
  and the user wants them packaged for a client.

## Slide types supported

| type         | Fields                                         |
| ------------ | ---------------------------------------------- |
| `title`      | `title`, `subtitle`                            |
| `section`    | `title`                                        |
| `bullets`    | `title`, `bullets[]`                           |
| `two_column` | `title`, `left_title`, `left_bullets[]`, `right_title`, `right_bullets[]` |
| `table`      | `title`, `headers[]`, `rows[][]`               |
| `image`      | `title`, `image_path`                          |

## Tests

```bash
pytest test_generate.py
```

The tests build a sample deck in-process and re-read it with `python-pptx`
to confirm each slide type renders.

## Guardrails

- **No client data in the spec files committed to this repo.** The YAML
  example above is placeholder.
- **Branding is neutral.** This skill does not ship a branded template.
  Projects that want custom branding should supply a `--template` path to a
  `.pptx` master; see the CLI help.
- **Accessibility.** Slide text is set as real text, not baked into images,
  so screen readers work.

## Extending

- Chart slides via `python-pptx` native charts (bar, line, pie) driven by
  data in the spec.
- Speaker notes via a `notes:` field on any slide type.
- Page numbers / footer logos via a template `.pptx`.
