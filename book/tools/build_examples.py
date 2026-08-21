#!/usr/bin/env python3
"""Run every book example through the real ``./Mathilda`` and emit its transcript.

For each ``book/examples/<path>.m`` this writes ``book/generated/<path>.tex``
containing the verified ``In[k]:= ... / Out[k]= ...`` transcript, which the book
includes verbatim via ``\\mtranscript{<path>}``.

Design contract (the reason this tool exists):
  * The ONLY source of an ``Out[]`` line is the binary.  Nothing here or in the
    ``.tex`` sources can hand-write an output, so no output can be guessed or
    silently drift away from what Mathilda actually prints.
  * Each ``.m`` file is ONE Mathilda session -- state carries across its lines,
    so an example may assign on one line and use the value on the next.  Session
    scope is therefore explicit and author-controlled: one file, one session.

Input format: one input expression per line; blank lines and lines beginning
with ``#`` are ignored.  A line whose result is empty (a ``;``-suppressed or
``Null`` setup line) prints its ``In[]`` with no ``Out[]`` line, matching the
convention used by the site's tutorials.

The NDJSON pipe driver is reused verbatim from ``site/verify_tutorial.py`` so the
book and the site stay byte-for-byte consistent.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]          # repo root
sys.path.insert(0, str(ROOT / "site"))
from verify_tutorial import run_session, MATHILDA    # noqa: E402  (reuse exact driver)

EXDIR = ROOT / "book" / "examples"
GENDIR = ROOT / "book" / "generated"


def render_pair(k, inp, out):
    """One In[k]:= / Out[k]= pair as plain text (Out omitted when empty/Null)."""
    lines = [f"In[{k}]:= {inp}"]
    # A ';'-suppressed or Null-valued setup line shows its In[] but no Out[],
    # matching the site tutorials' convention (generate.py drops Null too).
    if out.strip() not in ("", "Null"):
        lines.append(f"Out[{k}]= {out}")
    return "\n".join(lines) + "\n"


def render(inputs, outs):
    """The whole transcript (all pairs, blank line between) for \\mtranscript."""
    blocks = [render_pair(k, inp, out).rstrip()
              for k, (inp, out) in enumerate(zip(inputs, outs), 1)]
    return "\n\n".join(blocks) + "\n"


def read_inputs(mfile):
    return [l.rstrip("\n") for l in mfile.read_text().splitlines()
            if l.strip() and not l.lstrip().startswith("#")]


def main():
    if not MATHILDA.exists():
        sys.exit(f"error: {MATHILDA} not found -- run `make` in the repo root first.")
    mfiles = sorted(EXDIR.rglob("*.m"))
    if not mfiles:
        print(f"build_examples: no example files under {EXDIR.relative_to(ROOT)}")
        return
    total = 0
    for mf in mfiles:
        inputs = read_inputs(mf)
        outs = run_session(inputs)
        if len(outs) != len(inputs):
            print(f"build_examples: !! {mf.name}: {len(inputs)} inputs, "
                  f"{len(outs)} outputs", file=sys.stderr)
        rel = mf.relative_to(EXDIR)
        dst = (GENDIR / rel).with_suffix(".tex")
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(render(inputs, outs))
        # Per-pair snippets for \pair{...}: generated/<path>/<k>.tex holds In[k]/Out[k].
        pdir = (GENDIR / rel).with_suffix("")
        pdir.mkdir(parents=True, exist_ok=True)
        for k, (inp, out) in enumerate(zip(inputs, outs), 1):
            (pdir / f"{k}.tex").write_text(render_pair(k, inp, out))
        total += 1
        print(f"build_examples: {mf.relative_to(ROOT)} -> "
              f"{dst.relative_to(ROOT)} ({len(inputs)} inputs)")
    print(f"build_examples: generated {total} transcript(s).")


if __name__ == "__main__":
    main()
