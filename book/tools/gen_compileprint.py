#!/usr/bin/env python3
"""Generate CompilePrint bytecode listings for the book's \\compileprint boxes.

``CompilePrint[cf]`` disassembles a compiled function to the terminal via a raw
stdout write and returns ``Null`` -- so, unlike an ordinary ``In[]/Out[]`` result,
it is NOT carried on the NDJSON pipe that ``build_examples.py`` reads (that driver
skips every non-JSON line, exactly where the disassembly lives). This tool is the
capture mechanism for that output, in the same spirit as ``gen_figures.py`` for
plots: the book's verified-example promise (CONTEXT Sec. 3) is kept because the
listing is still produced by the real binary at build time, never by hand.

For each ``book/compileprint/<path>.m`` -- one file, one ``Compile[...]``
expression -- it runs ``CompilePrint[<that expression>]`` in its own Mathilda
process, captures the raw (non-JSON) stdout, and writes it verbatim to
``book/generated/compileprint/<path>.txt``, which ``\\compileprint{<path>}``
(book/mathilda.sty) renders in a titled box.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "site"))
from verify_tutorial import MATHILDA  # noqa: E402 (reuse the resolved binary path)

BOOK = ROOT / "book"
SRC = BOOK / "compileprint"
OUT = BOOK / "generated" / "compileprint"

# listings under inputenc utf8 rejects multibyte bytes; the disassembly is ASCII
# today, but sanitise defensively so a future opcode rendering can never abort the
# build. Mirrors gen_usage._ascii_safe.
_ASCII_MAP = {
    "—": "--", "–": "-", "−": "-", "…": "...",
    "→": "->", "←": "<-", "≤": "<=", "≥": ">=",
    "×": "x", "·": ".",
}


def _ascii_safe(s):
    for k, v in _ASCII_MAP.items():
        s = s.replace(k, v)
    return s.encode("ascii", "ignore").decode("ascii")


def read_expr(path):
    """One Compile[...] expression: non-comment, non-blank lines joined by space."""
    parts = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if s and not s.startswith("#"):
            parts.append(s)
    return " ".join(parts)


def disassemble(expr):
    """Run CompilePrint[<expr>] in its own process; return the raw disassembly
    (every non-JSON stdout line), stripped of trailing blank lines."""
    env = dict(os.environ)
    env["MATHILDA_NO_WINDOW"] = "1"
    proc = subprocess.Popen([str(MATHILDA)], stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                            text=True, bufsize=1, env=env)
    proc.stdin.write(json.dumps({"id": 1, "expr": f"CompilePrint[{expr}]"}) + "\n")
    proc.stdin.write('{"type":"quit"}\n')
    proc.stdin.flush()
    proc.stdin.close()
    raw = []
    for line in proc.stdout:
        s = line.rstrip("\n")
        st = s.strip()
        if st.startswith("{") and '"type"' in st:
            continue  # NDJSON control/result line -- not part of the disassembly
        raw.append(s)
    proc.wait(timeout=120)
    while raw and not raw[0].strip():
        raw.pop(0)
    while raw and not raw[-1].strip():
        raw.pop()
    return "\n".join(raw)


def main():
    if not MATHILDA.exists():
        sys.exit(f"error: {MATHILDA} not found -- run `make` in the repo root first.")
    srcs = sorted(SRC.rglob("*.m"))
    if not srcs:
        print("gen_compileprint: no book/compileprint/**/*.m sources found")
        return
    for src in srcs:
        rel = src.relative_to(SRC).with_suffix("")
        expr = read_expr(src)
        text = _ascii_safe(disassemble(expr))
        dst = OUT / f"{rel}.txt"
        dst.parent.mkdir(parents=True, exist_ok=True)
        if not text.strip():
            text = f"[CompilePrint produced no output for: {expr}]"
        dst.write_text(text + "\n")
        print(f"gen_compileprint: {rel} ({len(text.splitlines())} lines)")
    print(f"gen_compileprint: wrote {len(srcs)} listing(s) -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
