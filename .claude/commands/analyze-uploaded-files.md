---
description: Extract, summarize, or tabulate enterprise documents (.docx, .xlsx, .pdf, images, etc.) that the user has shared.
argument-hint: <file-path-or-directory>
allowed-tools: Read, Bash, Glob, Grep
---

The user has shared one or more files for you to analyze. The path(s) are
in `$ARGUMENTS`. If empty, list attached files or ask the user.

## Step 1 — Inventory

List the files and their types:

!`file $ARGUMENTS 2>/dev/null || ls -la $ARGUMENTS`

## Step 2 — Extract

For each file, use the [`analyze-enterprise-documents`](../skills/analyze-enterprise-documents/SKILL.md)
skill's dispatch table to pick the right extraction library:

- `.docx` → `python-docx`
- `.xlsx` / `.xls` → `openpyxl` or `pandas.read_excel`
- `.pdf` → `pypdf`, falling back to OCR for scans
- `.png` / `.jpg` → native multimodal read, or `pytesseract` for OCR-only
- `.csv` → `pandas.read_csv`
- `.pptx` → `python-pptx`

For large files, chunk — don't pour everything into one prompt.

## Step 3 — Reason

Based on what the user asked for, do one of:

- **Summarize** — the document's key points, structure, and any tables.
- **Tabulate** — pull numeric data into a structured format the user can
  copy.
- **Cross-reference** — compare claims across two or more files.
- **Answer specific questions** — e.g., "what does this PDF say about
  revenue in Q3?"

## Step 4 — Report with provenance

Every finding cites the source: file name, page / sheet / slide. If you
extracted a table, show its shape and source.

## Guardrails

See [`../shared/guidelines/documents-and-data.md`](../shared/guidelines/documents-and-data.md):

- **Do not send extracted content to external services** (translation APIs,
  pastebins, public LLMs beyond the current conversation) without asking.
- **Flag PII.** If you notice client identifiers, employee data, or
  financial figures, mention it before quoting at length.
- **Respect size limits** — chunk large PDFs / spreadsheets.

## Output shape

```
## Files analyzed
- report.pdf (14 pages, pypdf)
- figures.xlsx (3 sheets: Summary, Q3, Q4, openpyxl)

## Key findings
1. <finding, cited back to source>
2. ...

## Tables / data
<structured output if applicable>

## Flags
- PII detected in: <file / location>
- <any extraction quality issues>
```

## Related

- Skill: [`../skills/analyze-enterprise-documents/SKILL.md`](../skills/analyze-enterprise-documents/SKILL.md)
- Agent: [`../agents/enterprise-document-analyst.md`](../agents/enterprise-document-analyst.md)
  — dispatch this agent for large, multi-document analyses.
