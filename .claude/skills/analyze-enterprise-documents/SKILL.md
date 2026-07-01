---
name: analyze-enterprise-documents
description: Use when the user references or attaches Microsoft Office files (.doc, .docx, .xls, .xlsx, .ppt, .pptx), PDFs, images (.png, .jpg, .jpeg), or CSV/JSON/YAML from Teams or SharePoint and asks to extract, summarize, tabulate, or cross-reference their contents. Ships a runnable `extract.py` — invoke that directly rather than re-implementing extraction inline.
---

# analyze-enterprise-documents

**Ships runnable code.** This skill is not just a prompt — it includes
[`extract.py`](extract.py), a CLI + importable function that handles the
full matrix of common enterprise formats.

## Use it directly

```bash
# Install deps (once)
pip install -r requirements.txt

# CLI
python extract.py /path/to/report.pdf --max-pages 20
python extract.py /path/to/deck.pptx --format text
python extract.py /path/to/financials.xlsx > out.json

# Or import it
from extract import extract
e = extract("/path/to/report.pdf", max_pages=20)
print(e.type, e.metadata, len(e.content["pages"]))
```

Output is a structured `Extraction` dataclass with `type`, `source`,
`method`, `content`, `metadata`, and `flags` (e.g. `pdf_appears_scanned`,
`truncated_to_first_N_pages`).

## Supported formats

| Extension         | Library                       | Notes                                    |
| ----------------- | ----------------------------- | ---------------------------------------- |
| `.docx`           | `python-docx`                 | paragraphs + tables                      |
| `.xlsx` / `.xls`  | `pandas` + `openpyxl`         | multi-sheet; preserves types             |
| `.pdf`            | `pypdf`                       | flag set if pages appear scanned         |
| `.pptx`           | `python-pptx`                 | title + body + speaker notes per slide   |
| `.csv`            | `pandas`                      | all columns / rows                       |
| `.json` / `.yaml` | stdlib / PyYAML               | parsed data                              |
| `.txt` / `.md`    | `read_text`                   | raw text + metadata                      |
| `.png` / `.jpg`   | `pytesseract` (optional)      | OCR fallback; prefer multimodal read     |

## When an agent invokes this skill

1. Call `extract.extract(path)` instead of re-implementing.
2. Check `flags` — act on `pdf_appears_scanned_or_empty` by asking the user
   whether to OCR. Act on `truncated_to_first_N_pages` by chunking further.
3. Cite provenance in downstream output: include `result.source` and
   `result.method` in whatever the agent hands back to the user.

## Tests

```bash
pytest test_extract.py
```

The tests generate fixtures in-process (no committed binaries) to verify
each format handler. `reportlab` is imported on demand for the PDF test;
`pip install reportlab` if you want full coverage locally.

## Guardrails

See [`../../shared/guidelines/documents-and-data.md`](../../shared/guidelines/documents-and-data.md):

- **Do not** send extracted content to third-party services without approval.
- **Flag PII** on the way through — the skill is neutral; the calling agent
  decides how to redact.
- **Chunk large inputs** — pass `max_pages` for big PDFs; iterate over
  sheet subsets for wide spreadsheets.

## Extending

- Legacy `.doc`: add a `LibreOffice` fallback (`soffice --headless --convert-to docx`).
- Scanned PDFs: layer `pdf2image` + `pytesseract` under `_extract_pdf`.
- `.xlsb`: add `pyxlsb` to the xlsx dispatch branch.
- `.msg` (Outlook): add `extract-msg`.
