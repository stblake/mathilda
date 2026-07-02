---
name: generate-xlsx-deliverable
description: Use when the user asks to produce a formatted Excel workbook — scorecards, data tables, or multi-sheet reports — with consistent styling (headers, number formats, freeze panes, zebra striping). Ships a runnable `generate.py` plus a Python API.
---

# generate-xlsx-deliverable

Turns a spec (or a `dict[str, pandas.DataFrame]`) into a styled multi-sheet
`.xlsx` file. Takes the hand-formatting out of the loop so data analysis
lands in Excel without manual work.

## Use it directly

### From the CLI with a spec file

```bash
pip install -r requirements.txt
python generate.py spec.yaml --output deliverable.xlsx
```

Minimal spec:

```yaml
sheets:
  - name: Summary
    title: Q3 performance summary
    columns: [Metric, Q3 Actual, Q3 Plan, Delta]
    rows:
      - [Revenue ($M), 124, 118, 6]
      - [Costs ($M), 56, 58, -2]
      - [Margin (%), 55, 51, 4]
    number_formats:
      Q3 Actual: "#,##0"
      Q3 Plan: "#,##0"
      Delta: "+#,##0;-#,##0"
  - name: Detail
    title: Detail by segment
    columns: [Segment, Customers, ARR ($K)]
    rows:
      - [Enterprise, 42, 8400]
      - [Mid-market, 118, 6200]
      - [SMB, 430, 2100]
```

### From Python with dataframes

```python
from generate import write_workbook
import pandas as pd

sheets = {
    "Summary": pd.DataFrame({
        "Metric": ["Revenue", "Costs", "Margin"],
        "Q3": [124, 56, 55],
        "Plan": [118, 58, 51],
    }),
}
write_workbook(sheets, "deliverable.xlsx", title_row=True)
```

## Styling applied

- **Header row:** bold, filled, white text.
- **Freeze panes** just below the header.
- **Zebra striping** on data rows.
- **Column auto-width** capped to a sane maximum.
- **Number formats** per-column when provided.
- **Title row** above the header when a sheet's `title` is set.

Recipients recognize this shape; no further hand-formatting required for most
status / scorecard / data-table outputs.

## When to invoke

- User asks to "make an Excel", "export to xlsx", or "produce a summary workbook".
- Another skill produces `dict[str, DataFrame]` (e.g. analysis on top of
  [`analyze-enterprise-documents`](../analyze-enterprise-documents/SKILL.md))
  and the user wants it packaged.

## Tests

```bash
pytest test_generate.py
```

Writes sample workbooks to a temp dir and re-reads them with openpyxl
to confirm styling and content.

## Guardrails

- **No client data in committed sample specs.** The YAML above is
  placeholder.
- **Numbers preserved exactly.** No rounding on write — use
  `number_formats` for display.
- **No formulas generated** — this skill writes values only. If the client
  wants formulas, the producing code should pass them as strings starting
  with `=`, and the caller accepts they won't be recomputed by openpyxl.

## Extending

- Conditional formatting (data bars, color scales) by column.
- Native Excel charts (line, bar) via openpyxl's chart API.
- Template `.xlsx` with client branding applied as a starting workbook.
