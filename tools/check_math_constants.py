#!/usr/bin/env python3
"""Guard against unguarded POSIX <math.h> constants (M_PI, M_E, ...).

`M_PI`, `M_E` and friends are POSIX, *not* C99. glibc hides them under
`-std=c99`, so a file that uses one without a fallback fails to compile on
Linux. Darwin exposes them implicitly, which masks the bug locally and lets it
reach users as a build failure (see issue #36).

The convention -- documented in CLAUDE.md and SPEC.md Sec. 10, and followed by
src/trig.c and src/numeric.c -- is a per-file guarded fallback placed right
after `#include <math.h>`:

    #include <math.h>
    #ifndef M_PI
    #define M_PI 3.14159265358979323846
    #endif

This script enforces that convention across the tree. Run it via
`make check-c99`; it exits nonzero (and names every offending file) when a
constant is used without a matching guard in the same translation unit.

Guards are checked per *file*, not per translation unit reachable through
includes, because that is the convention: no project header defines these.
"""

import os
import re
import subprocess
import sys

# The full set of POSIX <math.h> constants. All carry the same portability
# hazard, so all are checked -- not just the ones that have bitten us so far.
CONSTANTS = (
    "M_E", "M_LOG2E", "M_LOG10E", "M_LN2", "M_LN10",
    "M_PI", "M_PI_2", "M_PI_4", "M_1_PI", "M_2_PI",
    "M_2_SQRTPI", "M_SQRT2", "M_SQRT1_2",
)

# Reference values, so the error message can be pasted straight into the fix.
VALUES = {
    "M_E":        "2.71828182845904523536",
    "M_LOG2E":    "1.44269504088896340736",
    "M_LOG10E":   "0.434294481903251827651",
    "M_LN2":      "0.693147180559945309417",
    "M_LN10":     "2.30258509299404568402",
    "M_PI":       "3.14159265358979323846",
    "M_PI_2":     "1.57079632679489661923",
    "M_PI_4":     "0.785398163397448309616",
    "M_1_PI":     "0.318309886183790671538",
    "M_2_PI":     "0.636619772367581343076",
    "M_2_SQRTPI": "1.12837916709551257390",
    "M_SQRT2":    "1.41421356237309504880",
    "M_SQRT1_2":  "0.707106781186547524401",
}

USE_RE = re.compile(r"\b(" + "|".join(CONSTANTS) + r")\b")
GUARD_RE = re.compile(r"^\s*#\s*(?:ifndef|define)\s+(M_\w+)")


def strip_comments_and_strings(text):
    """Blank out comments and string/char literals.

    A constant named only in prose ("M_PI is POSIX, not C99") is not a usage,
    and a checker that cannot tell the difference gets ignored. Replacing with
    spaces rather than deleting keeps line numbers intact for reporting.
    """
    out = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        two = text[i:i + 2]
        if two == "/*":
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            out.append("".join(ch if ch == "\n" else " " for ch in text[i:j]))
            i = j
        elif two == "//":
            j = text.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i))
            i = j
        elif c in "\"'":
            quote, j = c, i + 1
            while j < n and text[j] != quote:
                j += 2 if text[j] == "\\" else 1
            j = min(j + 1, n)
            out.append("".join(ch if ch == "\n" else " " for ch in text[i:j]))
            i = j
        else:
            out.append(c)
            i += 1
    return "".join(out)


def source_files(repo_root):
    """Git-tracked C sources under src/ and tests/, excluding vendored code."""
    try:
        tracked = subprocess.check_output(
            ["git", "-C", repo_root, "ls-files", "src", "tests"],
            text=True,
        ).split("\n")
    except (subprocess.CalledProcessError, OSError) as exc:
        sys.stderr.write("check_math_constants: cannot list git files: %s\n" % exc)
        sys.exit(2)

    ext = (".c", ".h", ".inc")
    return [
        f for f in tracked
        if f.endswith(ext) and not f.startswith("src/external/")
    ]


def scan(repo_root, path):
    """Return the sorted constants used in `path` that it does not guard."""
    try:
        with open(os.path.join(repo_root, path), encoding="utf-8",
                  errors="replace") as fh:
            raw = fh.read()
    except OSError:
        return []

    used, guarded = set(), set()
    for line in strip_comments_and_strings(raw).splitlines():
        guard = GUARD_RE.match(line)
        if guard:
            # `#ifndef M_PI` / `#define M_PI ...` -- the guard itself, not a use.
            guarded.add(guard.group(1))
            continue
        used.update(USE_RE.findall(line))
    return sorted(used - guarded)


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    offenders = []
    for path in source_files(repo_root):
        missing = scan(repo_root, path)
        if missing:
            offenders.append((path, missing))

    if not offenders:
        return 0

    sys.stderr.write(
        "error: POSIX <math.h> constants used without a C99 fallback.\n"
        "       glibc hides these under -std=c99; the build fails on Linux\n"
        "       even though macOS compiles it. Add, after #include <math.h>:\n\n"
    )
    for path, missing in offenders:
        sys.stderr.write("  %s\n" % path)
        for name in missing:
            sys.stderr.write("      #ifndef %s\n      #define %s %s\n      #endif\n"
                             % (name, name, VALUES[name]))
    sys.stderr.write("\n%d file(s) affected.\n" % len(offenders))
    return 1


if __name__ == "__main__":
    sys.exit(main())
