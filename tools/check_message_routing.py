#!/usr/bin/env python3
"""
check_message_routing.py -- does every user-facing message support Quiet[] and
Check[]?

THE FAILURE MODE THIS EXISTS FOR. Mathilda's Quiet[] (suppression) and Check[]
(detection) are not a single choke-point: they are a two-step convention every
emission site must uphold by hand --

    (1) mth_msg_note_fired()          so an enclosing Check[] sees it, and
    (2) if (mth_msg_suppressed()) ..  so Quiet[] silences it.

A diagnostic printed with a RAW `fprintf(stderr, "Head::tag: ...")` upholds
neither, so Quiet[] does not suppress it AND Check[] does not detect it. That is
a WRONG answer for Check[Head[bad], fallback] (it returns Head[bad], not the
fallback) and a leak past Quiet[] -- not merely cosmetic. It shipped that way at
~340 sites before this gate: MatrixPower escaped both, and Power::infy /
Infinity::indet were invisible to Check[].

THE FIX is a single funnel, src/message.c's mth_message() family, that does both
steps once. A migrated site reads `mth_message("Power", "infy", "...")` -- head
and tag are SEPARATE arguments, so no "::"-bearing literal ever reaches an
`fprintf(stderr, ...)`. This gate keys off the raw stderr write, so a routed
site simply vanishes from its scan; only a bypassing site is named.

WHY IT READS THE SOURCE. Dispatch lives in the source, and a message that
bypasses the funnel does so silently -- nothing fails, Check just quietly misses
it. So this diffs the raw `Head::tag` stderr writes in src/ against a checked-in
EXEMPT (permanent, non-user-facing) + BASELINE (the shrinking migration
backlog). It RATCHETS: it fails on a NEW bypassing site and on a BASELINE entry
that is no longer detected (a migrated site whose entry must now be deleted), so
the backlog cannot quietly grow and cannot quietly rot.

A site may legitimately not route through the funnel: the funnel body itself,
parser syntax errors (they run before any Check[] frame exists), the REPL
file:line reporter, and OOM aborts. Those go in EXEMPT WITH A REASON. Everything
else is BASELINE until migrated, then removed.

Usage:  python3 tools/check_message_routing.py     (exit 1 on a finding)
        make check-messages
"""

import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")

# ---------------------------------------------------------------------------
# Permanently exempt. A key is either a whole file (relpath) or a single
# `relpath:Head::tag` site. Every entry needs a reason -- the reason is the
# point of the list, so a reader can tell "considered and rejected" from "never
# noticed". These are NOT user-facing math messages, or they ARE the funnel.
# ---------------------------------------------------------------------------
EXEMPT = {
    "src/message.c": "the funnel itself -- mth_message_v's `fprintf(stderr, "
    '"%s::%s: ", head, tag)` is the one legitimate raw Head::tag stderr write, '
    "and builtin_message prints a user Message[]'s already-resolved template "
    "here. This IS the choke-point every other site routes through.",
}

# ---------------------------------------------------------------------------
# The migration backlog: raw Head::tag stderr sites not yet routed through the
# funnel. Keyed `relpath:Head::tag` (a `%s` head/tag is a shared helper that
# fans out to many heads). RATCHET: shrinks as each module migrates; a stale
# entry (no longer detected) fails the gate so it is deleted when its site is
# routed. Target end-state: empty (assert-empty), kept as a visible debt marker.
# Seeded from the initial audit; see tasks/todo.md / the plan.
# ---------------------------------------------------------------------------
BASELINE = {
}


def read(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        return f.read()


def all_sources():
    for dirpath, _dirs, files in os.walk(SRC):
        if os.sep + "external" in dirpath:
            continue
        for fn in files:
            if fn.endswith((".c", ".h")):
                yield os.path.join(dirpath, fn)


def blank_comments(text):
    """Replace /* ... */ and // ... comments with equal-length spaces, keeping
    newlines, so (a) reported line numbers stay accurate and (b) `Head::tag`
    prose inside a comment or docstring is never mistaken for an emission site.
    A hand scanner (not a regex) so a `//` inside a string or a `"` inside a
    comment cannot desync the state -- string and char literals are left INTACT
    because the format argument we classify is a string literal."""
    out = []
    i, n = 0, len(text)
    state = "code"  # code | line | block | str | chr
    while i < n:
        c = text[i]
        nxt = text[i + 1] if i + 1 < n else ""
        if state == "code":
            if c == "/" and nxt == "/":
                state = "line"; out.append("  "); i += 2; continue
            if c == "/" and nxt == "*":
                state = "block"; out.append("  "); i += 2; continue
            if c == '"':
                state = "str"; out.append(c); i += 1; continue
            if c == "'":
                state = "chr"; out.append(c); i += 1; continue
            out.append(c); i += 1; continue
        if state == "line":
            if c == "\n":
                state = "code"; out.append("\n"); i += 1; continue
            out.append(" "); i += 1; continue
        if state == "block":
            if c == "*" and nxt == "/":
                state = "code"; out.append("  "); i += 2; continue
            out.append("\n" if c == "\n" else " "); i += 1; continue
        if state in ("str", "chr"):
            if c == "\\":
                out.append(c)
                if nxt:
                    out.append(nxt)
                i += 2; continue
            if (state == "str" and c == '"') or (state == "chr" and c == "'"):
                state = "code"; out.append(c); i += 1; continue
            out.append(c); i += 1; continue
    return "".join(out)


# A stderr write: fprintf/gmp_fprintf(stderr, <fmt>, ...) or fputs("...", stderr).
STDERR_FPRINTF = re.compile(r"\b(?:fprintf|gmp_fprintf)\s*\(\s*stderr\s*,", re.S)
STDERR_FPUTS = re.compile(
    r'\bfputs\s*\(\s*("(?:\\.|[^"\\])*")\s*,\s*stderr\s*\)', re.S
)
# One-or-more adjacent string literals immediately after the comma (C
# concatenates them); this is the format argument.
ADJ_LITERALS = re.compile(r'\s*((?:"(?:\\.|[^"\\])*"\s*)+)')
ONE_LITERAL = re.compile(r'"((?:\\.|[^"\\])*)"')
# A Wolfram-style diagnostic: Head::tag or $Head::tag or %s::(tag|%s), anchored
# at the start of the (concatenated) format string.
DIAGNOSTIC = re.compile(r"^(%s|[A-Za-z$][A-Za-z0-9$]*)::(%s|[a-z][A-Za-z0-9]*)")


def concat_literals(blob):
    """Join the contents of adjacent "..." string literals into one string."""
    return "".join(m.group(1) for m in ONE_LITERAL.finditer(blob))


def line_of(text, pos):
    return text.count("\n", 0, pos) + 1


def detect(text):
    """Yield (Head::tag key-suffix, line) for every raw diagnostic stderr write.
    Operates on comment-blanked text so line numbers match the original."""
    for m in STDERR_FPRINTF.finditer(text):
        lit = ADJ_LITERALS.match(text, m.end())
        if not lit:
            continue
        content = concat_literals(lit.group(1))
        d = DIAGNOSTIC.match(content)
        if not d:
            continue
        head = d.group(1)
        tag = "<dynamic>" if d.group(2) == "%s" else d.group(2)
        yield f"{head}::{tag}", line_of(text, m.start())
    for m in STDERR_FPUTS.finditer(text):
        content = concat_literals(m.group(1))
        d = DIAGNOSTIC.match(content)
        if not d:
            continue
        head = d.group(1)
        tag = "<dynamic>" if d.group(2) == "%s" else d.group(2)
        yield f"{head}::{tag}", line_of(text, m.start())


def main():
    detected = {}  # key "relpath:Head::tag" -> sorted list of line numbers
    for path in all_sources():
        rel = os.path.relpath(path, ROOT)
        text = blank_comments(read(path))
        for suffix, line in detect(text):
            detected.setdefault(f"{rel}:{suffix}", []).append(line)

    def exempt(key):
        # A key is "relpath:Head::tag"; relpaths never contain a colon, so the
        # first ":" splits the file off cleanly. Match a whole-file exemption
        # (relpath) or a single-site one (the full key).
        relpath = key.split(":", 1)[0]
        return relpath in EXEMPT or key in EXEMPT

    problems = sorted(
        k for k in detected if not exempt(k) and k not in BASELINE
    )
    stale = sorted(k for k in BASELINE if k not in detected)

    total_sites = sum(len(v) for v in detected.values())
    print(
        f"scanned src/ (excl external): {len(detected)} distinct raw Head::tag "
        f"stderr keys over {total_sites} sites"
    )
    print(f"EXEMPT={len(EXEMPT)}  BASELINE={len(BASELINE)}")

    failed = False
    if problems:
        failed = True
        print("\nmessage-routing audit FAILED -- new raw diagnostic(s):\n", file=sys.stderr)
        for k in problems:
            lines = ",".join(str(x) for x in sorted(detected[k]))
            print(
                f"  {k}  (lines {lines})\n"
                f"    A raw fprintf(stderr, \"...\") that bypasses Quiet[]/Check[].\n"
                f"    Route it through src/message.c's mth_message()/mth_message_gated(),\n"
                f"    or add it to EXEMPT/BASELINE in this script with a reason.\n",
                file=sys.stderr,
            )

    if stale:
        failed = True
        print(
            "\nmessage-routing audit FAILED -- stale BASELINE entrie(s) "
            "(migrated? then delete them):\n",
            file=sys.stderr,
        )
        for k in stale:
            print(f"  {k}", file=sys.stderr)
        print("", file=sys.stderr)

    if failed:
        print(
            f"{len(problems)} new, {len(stale)} stale. "
            "See docs/design/message_routing.md.",
            file=sys.stderr,
        )
        return 1

    if BASELINE:
        print(
            f"OK: no new raw diagnostics. {len(BASELINE)} still on the migration "
            "backlog (BASELINE)."
        )
    else:
        print(
            "OK: every user-facing message routes through the Quiet/Check funnel "
            "(BASELINE empty)."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
