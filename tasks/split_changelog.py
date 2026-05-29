#!/usr/bin/env python3
"""Split monthly Mathilda changelogs into weekly files keyed by Monday of the ISO week.

Reads:
    docs/spec/changelog/2026-04.md
    docs/spec/changelog/2026-05.md

Writes:
    docs/spec/changelog/<YYYY-MM-DD>.md   (one per Monday-week)

Rules:
- H2 sections beginning with `## ` are the splitting unit.
- The section's date is parsed from the trailing `(YYYY-MM-DD)` in the heading
  line. Sections without a trailing date inherit the most recently seen date
  in their source file (whichever direction — both forward and reverse
  chronological orderings appear).
- Within a weekly file, sections are emitted newest-first.
"""

from __future__ import annotations

import datetime
import os
import re
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
CHANGELOG_DIR = ROOT / "docs" / "spec" / "changelog"
SOURCES = [CHANGELOG_DIR / "2026-04.md", CHANGELOG_DIR / "2026-05.md"]

HEADING_DATE_RE = re.compile(r"\((\d{4})-(\d{2})-(\d{2})\)\s*$")
ANY_DATE_RE = re.compile(r"\b(\d{4})-(\d{2})-(\d{2})\b")


def week_start(d: datetime.date) -> datetime.date:
    """Return the Monday of d's ISO week."""
    return d - datetime.timedelta(days=d.weekday())


def parse_sections(path: Path) -> list[tuple[datetime.date, str]]:
    """Return [(section_date, section_text_including_heading), ...] in source order."""
    text = path.read_text()
    # Split on '## ' at line start; first chunk is preamble (drop).
    parts = re.split(r"(?m)^## ", text)
    if not parts:
        return []
    out: list[tuple[datetime.date | None, str]] = []
    for chunk in parts[1:]:
        full = "## " + chunk
        first_line = chunk.split("\n", 1)[0]
        m = HEADING_DATE_RE.search(first_line)
        d: datetime.date | None
        if m:
            d = datetime.date(int(m.group(1)), int(m.group(2)), int(m.group(3)))
        else:
            d = None
        out.append((d, full))

    # Forward fill undated sections from the previous dated entry; if the file
    # starts with an undated section, backfill from the first dated one.
    last: datetime.date | None = None
    filled: list[tuple[datetime.date, str]] = []
    for d, body in out:
        if d is None:
            if last is None:
                # Look ahead for the first dated entry.
                for d2, _ in out:
                    if d2 is not None:
                        last = d2
                        break
                if last is None:
                    raise RuntimeError(f"No dates found anywhere in {path}")
            filled.append((last, body))
        else:
            last = d
            filled.append((d, body))
    return filled


def main() -> int:
    by_week: dict[datetime.date, list[tuple[datetime.date, str, int]]] = defaultdict(list)
    seq = 0
    total_input = 0
    for src in SOURCES:
        sections = parse_sections(src)
        total_input += len(sections)
        for d, body in sections:
            seq += 1
            by_week[week_start(d)].append((d, body, seq))

    total_output = 0
    written: list[tuple[datetime.date, Path, int]] = []
    for monday, items in sorted(by_week.items()):
        # Newest first within the week. Stable secondary sort by original
        # sequence so same-date sections keep their relative order.
        items.sort(key=lambda t: (t[0], -t[2]), reverse=True)
        sunday = monday + datetime.timedelta(days=6)
        out_path = CHANGELOG_DIR / f"{monday.isoformat()}.md"
        header = (
            f"# Changelog: week of {monday.isoformat()} "
            f"(Mon) – {sunday.isoformat()} (Sun)\n\n"
            f"Feature additions and fixes recorded during this week.\n\n"
        )
        body = "".join(b for _, b, _ in items)
        if not body.endswith("\n"):
            body += "\n"
        out_path.write_text(header + body)
        total_output += len(items)
        written.append((monday, out_path, len(items)))

    print(f"Input H2 sections : {total_input}")
    print(f"Output H2 sections: {total_output}")
    print()
    print("Per-week output:")
    for monday, path, n in written:
        print(f"  {monday.isoformat()}.md  sections={n}  path={path.relative_to(ROOT)}")

    if total_input != total_output:
        print("MISMATCH — investigate", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
