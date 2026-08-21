#!/usr/bin/env python3
"""Verify that every ``\\B{Name}`` in the book points at a real reference page.

Scans the book's ``.tex`` sources for ``\\B{...}`` uses and checks each name
against ``site/docs/assets/builtins.json``.  Exits non-zero (and lists the
offenders) if any ``\\B{}`` names a builtin with no reference page -- catching
typos and links to symbols that were never documented.

Usage: python3 book/tools/check_links.py
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
JSON = ROOT / "site" / "docs" / "assets" / "builtins.json"
BOOK = ROOT / "book"

# \B{...} -- capture the braced argument (names never contain a nested brace).
B_RE = re.compile(r"\\B\{([^{}]*)\}")


def known_names():
    data = json.load(open(JSON))
    # Store both the raw name and the LaTeX-escaped ($ -> \$) form, since authors
    # write $-symbols as \B{\$Version} to keep the dollar out of math mode.
    names = set()
    for e in data:
        names.add(e["name"])
        names.add(e["name"].replace("$", r"\$"))
    return names


def scan_uses():
    uses = {}                                          # name -> [files]
    for tex in sorted(BOOK.rglob("*.tex")):
        if "generated" in tex.parts:
            continue
        for m in B_RE.finditer(tex.read_text()):
            uses.setdefault(m.group(1), []).append(tex.relative_to(ROOT))
    return uses


def main():
    if not JSON.exists():
        sys.exit(f"error: {JSON} not found.")
    known = known_names()
    uses = scan_uses()
    unlinked = {n: f for n, f in uses.items() if n not in known}
    total = sum(len(v) for v in uses.values())
    print(f"check_links: {total} \\B{{}} use(s), {len(uses)} distinct name(s).")
    if unlinked:
        print(f"check_links: FAIL -- {len(unlinked)} name(s) with no reference page:")
        for n, files in sorted(unlinked.items()):
            where = ", ".join(str(f) for f in dict.fromkeys(files))
            print(f"  \\B{{{n}}}  (in {where})")
        return 1
    print("check_links: OK -- every \\B{} resolves to a reference page.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
