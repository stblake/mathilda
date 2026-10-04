#!/usr/bin/env python3
"""Backstop gate: every example on a published doc page still matches the binary.

`make docs` verifies examples at *generation* time, so a freshly regenerated
tree is correct by construction. This gate guards the other direction: a
committed tree that has gone **stale** because the binary's behaviour changed
and nobody re-ran `make docs`. It re-extracts every ``In[]:=``/``Out[]=`` pair
from the generated pages under ``site/docs/documentation/`` and re-runs it
through the current ``./Mathilda``.

Each page is replayed in ONE session in document order (as a reader runs it top
to bottom, and as `generate.py` groups its examples), so a setup line that binds
a variable is visible to the later examples that use it. Inputs that carry no
``Out[]`` (``;``-suppressed / ``Null`` setup lines) are still fed, to keep that
state, but are not compared.

Exit status is nonzero if any non-exempt page has a mismatch.

Usage:
  python3 site/verify_docs_examples.py                  # whole tree
  python3 site/verify_docs_examples.py special-functions # one category
  python3 site/verify_docs_examples.py special-functions/EllipticF
"""
import concurrent.futures
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import generate  # noqa: E402  (reuse run_session + the correct binary path)

DOC = generate.DOC_OUT

IN_RE = re.compile(r"^In\[\d+\]:=\s?(.*)$")
OUT_RE = re.compile(r"^Out\[\d+\]=\s?(.*)$")
FENCE_RE = re.compile(r"```mathematica\n(.*?)```", re.S)
HEADING_RE = re.compile(r"^#{1,6} ", re.M)

# Pages this fast backstop cannot reproduce by replaying the rendered page.
# The authoritative verifier is `make docs`, which runs each example in the
# generator's exact session grouping and spec order; this gate re-extracts the
# *rendered* page, which (a) loses that grouping and (b) can reorder examples
# across subsections, so a handful of pages genuinely cannot be checked this way.
# Each entry is "<category>/<Name>" with a one-line reason. Keep it short.
EXEMPT = {
    # Nondeterministic output — unverifiable by any replay:
    "expression-information/MemoryInUse",  # live heap size varies run to run
    "graphics/Histogram",                  # RandomReal[] sample data
    "time-and-date/Timing",                # wall-clock time varies
    "time-and-date/RepeatedTiming",        # wall-clock time varies
    "time-and-date/AbsoluteTime",          # wall-clock seconds; sci-notation rounding drifts
    # Rendered order/grouping != execution order (generator runs spec blocks in
    # one session and in spec order; render_examples splits/reorders them):
    "linear-algebra/Inverse",       # Options example reuses a matrix set in Basic
    "linear-algebra/MatrixRank",    # Options examples depend on `m` set in Basic
    "machine-learning/LearnDistribution",  # Basic uses `d` defined later in Options
    "machine-learning/DimensionReduce",    # Options example reuses `d` set in Basic
}


def extract_pairs(text):
    """Return [(input, expected_or_None), ...] across all fences, in order.

    ``expected`` is None for a setup line (an input with no ``Out[]``); such a
    line is still fed to the session but not compared. A multi-line ``Out[]``
    (a matrix, a series) is captured until the next In/Out line, a blank line,
    or the end of the fence."""
    pairs = []
    for fence in FENCE_RE.findall(text):
        lines = fence.split("\n")
        i = 0
        while i < len(lines):
            m = IN_RE.match(lines[i])
            if not m:
                i += 1
                continue
            expr = m.group(1).strip()
            i += 1
            if i < len(lines) and OUT_RE.match(lines[i]):
                out_lines = [OUT_RE.match(lines[i]).group(1)]
                i += 1
                while i < len(lines):
                    ln = lines[i]
                    if ln.strip() == "" or IN_RE.match(ln) or OUT_RE.match(ln):
                        break
                    out_lines.append(ln)
                    i += 1
                pairs.append((expr, "\n".join(out_lines).strip()))
            elif expr:
                pairs.append((expr, None))
    return pairs


def split_sections(text):
    """Split a page at every heading line, so each ``###`` subsection is its own
    chunk. `generate.py` runs each example group (which renders as one
    subsection) through a SEPARATE binary session, so a symbol a Basic-examples
    line assigns is invisible to an Applications line below it. Replaying the
    whole page as one session would leak that state and report false mismatches;
    replaying each heading chunk as its own session mirrors the generator."""
    idxs = [m.start() for m in HEADING_RE.finditer(text)]
    if not idxs:
        return [text]
    bounds = idxs + [len(text)]
    parts = [text[:idxs[0]]] if idxs[0] > 0 else []
    parts += [text[idxs[i]:bounds[i + 1]] for i in range(len(idxs))]
    return parts


def check_page(page):
    """Return (rel, [mismatch, ...]). A mismatch is (input, expected, actual).

    Each heading-delimited section is replayed as its own session (see
    split_sections); within a section, inputs run in document order so a setup
    line binds for the lines that follow it."""
    rel = f"{page.parent.name}/{page.stem}"
    text = page.read_text(errors="replace")
    mismatches = []
    for section in split_sections(text):
        pairs = extract_pairs(section)
        inputs = [p[0] for p in pairs if p[0]]
        if not inputs:
            continue
        result = generate.run_session(inputs)
        if result is None:
            mismatches.append(("<session>", "", "TIMEOUT"))
            continue
        outs, _ = result
        j = 0
        for expr, expected in pairs:
            if not expr:
                continue
            actual = (outs[j] if j < len(outs) else "").strip()
            j += 1
            if expected is None:
                continue
            if actual != expected.strip():
                mismatches.append((expr, expected.strip(), actual))
    return rel, mismatches


def main():
    args = sys.argv[1:]
    if args:
        pages = []
        for a in args:
            p = DOC / (a + ".md") if "/" in a else None
            if p and p.exists():
                pages.append(p)
            else:
                pages.extend(sorted((DOC / a).glob("*.md")))
        pages = [p for p in pages if p.name != "index.md"]
    else:
        pages = [p for p in sorted(DOC.glob("*/*.md")) if p.name != "index.md"]

    # Examples that do real filesystem I/O against shared /tmp paths can race
    # when run concurrently (one page's write truncates a file another page is
    # reading). This is not only the file-io category — other-advanced carries
    # the Read-family (Word, EndOfFile, Character, ...) which also round-trips
    # through /tmp. Detect I/O pages by content and run them serially (each is a
    # self-contained create-then-read, so serial execution is race-free);
    # everything else goes through the thread pool.
    io_markers = re.compile(
        r"Open(Write|Append|Read)|WriteString|\bWrite\[|\bExport\[|\bPut\[|"
        r"DeleteFile|/tmp/")
    def is_io(p):
        return (p.parent.name == "file-io"
                or bool(io_markers.search(p.read_text(errors="replace"))))
    io_pages = [p for p in pages if is_io(p)]
    par_pages = [p for p in pages if p not in io_pages]

    workers = min(os.cpu_count() or 4, 8)
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        results.extend(ex.map(check_page, par_pages))
    results.extend(check_page(p) for p in io_pages)   # serial, race-free

    bad, n_pages, n_ex = [], 0, 0
    for rel, mism in results:
        n_pages += 1
        if rel in EXEMPT:
            continue
        if mism:
            bad.append((rel, mism))

    for rel, mism in sorted(bad):
        for expr, expected, actual in mism:
            n_ex += 1
            print(f"MISMATCH  {rel}")
            print(f"    In   : {expr}")
            print(f"    page : {expected!r}")
            print(f"    build: {actual!r}")
    if bad:
        print(f"\nFAIL: {n_ex} mismatch(es) on {len(bad)} page(s) "
              f"of {n_pages} checked. Run 'make docs' and commit site/docs/**.")
        sys.exit(1)
    print(f"OK: {n_pages} pages checked, every example matches the current build.")


if __name__ == "__main__":
    main()
