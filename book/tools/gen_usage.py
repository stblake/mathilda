#!/usr/bin/env python3
"""Generate docstring ("usage") cards for builtins the book documents with \\usagebox.

Scans the book's ``.tex`` for ``\\usagebox{Name}``, asks the running ``./Mathilda``
for each name's usage string (``?Name``) and its ``Attributes``, and writes, under
``book/generated/usage/``:

  <Name>.txt   -- the raw usage text, exactly what ``?Name`` prints (verbatim body)
  <Name>.attr  -- the attribute list, comma-separated, no braces

so ``\\usagebox`` (book/mathilda.sty) can render a titled card. Nothing is
hand-written: the card body is the binary's own docstring.
"""
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "site"))
from verify_tutorial import MATHILDA  # noqa: E402 (reuse the resolved binary path)

BOOK = ROOT / "book"
OUT = ROOT / "book" / "generated" / "usage"
USAGE_RE = re.compile(r"\\usagebox\{([^{}]*)\}")


def collect_names():
    names = []
    for tex in sorted(BOOK.rglob("*.tex")):
        if "generated" in tex.parts:
            continue
        for m in USAGE_RE.finditer(tex.read_text()):
            if m.group(1) not in names:
                names.append(m.group(1))
    return names


def query(names):
    """One session: for each name send ?Name (usage) and Attributes[Name]."""
    lines, idmap, k = [], {}, 1
    for n in names:
        lines.append(json.dumps({"id": k, "expr": f"?{n}"})); idmap[k] = ("usage", n); k += 1
        lines.append(json.dumps({"id": k, "expr": f"Attributes[{n}]"})); idmap[k] = ("attr", n); k += 1
    inp = "\n".join(lines) + '\n{"type":"quit"}\n'
    env = dict(os.environ); env["MATHILDA_NO_WINDOW"] = "1"
    r = subprocess.run([str(MATHILDA)], input=inp, capture_output=True, text=True,
                       timeout=180, env=env)
    usage, attrs = {}, {}
    for line in r.stdout.splitlines():
        line = line.strip()
        if not (line.startswith("{") and '"type"' in line):
            continue
        try:
            m = json.loads(line)
        except ValueError:
            continue
        if "id" not in m:
            continue
        kind, name = idmap.get(m["id"], (None, None))
        if kind == "usage" and m.get("type") == "usage":
            usage[name] = m.get("payload", "")
        elif kind == "attr" and m.get("type") == "expr":
            attrs[name] = m.get("payload", "")
    return usage, attrs


def main():
    names = collect_names()
    if not names:
        print("gen_usage: no \\usagebox{...} uses found")
        return
    if not MATHILDA.exists():
        sys.exit(f"error: {MATHILDA} not found -- run `make` in the repo root first.")
    OUT.mkdir(parents=True, exist_ok=True)
    usage, attrs = query(names)
    for n in names:
        u = usage.get(n, "").rstrip("\n")
        (OUT / f"{n}.txt").write_text((u + "\n") if u else f"(no usage string for {n})\n")
        a = re.sub(r"^\s*\{|\}\s*$", "", attrs.get(n, "").strip())
        (OUT / f"{n}.attr").write_text(a)
        print(f"gen_usage: {n} ({len(u)} chars; attributes: {a})")
    print(f"gen_usage: wrote {len(names)} usage card(s) -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
