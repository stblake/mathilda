#!/usr/bin/env python3
"""check_bronstein_test_coverage.py -- every routine of the Bronstein campaign
has a direct unit test (BRONSTEIN_PLAN.md "Testing policy").

The routines are:
  * every `name::usage` in src/internal/bronstein/*.m;
  * the Trager routines this campaign added or generalised (TRAGER_NAMES);
  * the C builtins it added (C_NAMES).

A routine is covered when its name occurs as a whole word in one of the test
suites (TEST_GLOBS), outside a comment. The gate is assert-empty: a routine
can be left out only by an EXEMPT entry, which needs a one-line reason.
Stale EXEMPT entries (the name is now tested, or no longer exists) fail too,
so the list cannot rot.

Usage:  python3 tools/check_bronstein_test_coverage.py      (exit 1 on a gap)
"""
import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

TEST_GLOBS = [
    "tests/test_bronstein_*.c",
    "tests/test_integrate_cascade_trager_bronstein.c",
    "tests/test_trager_*.c",
    "tests/test_risch_rde_tower.c",
]

# Trager routines added or generalised for the rational-function field (Field.m)
TRAGER_NAMES = [
    "fieldFn", "kcFn", "rcFn", "fnInvertT", "fnCollect", "fieldExtQ",
    "fieldScalarSymbols", "fieldScalarQ", "fieldBase",
]

# C builtins and public symbols added by the campaign
C_NAMES = [
    "Integrate`Bronstein", "Integrate`CascadeGateQ", "Risch`RischDECertified",
    "Risch`RischDECertifiedReason", "$TragerCascadeBudget", "$BronsteinCascadeBudget",
]

# name -> reason (one line).  Keep this short.
EXEMPT = {
}


def strip_c_comments(src):
    src = re.sub(r"/\*.*?\*/", " ", src, flags=re.S)
    return re.sub(r"//[^\n]*", " ", src)


def main():
    names = {}
    for path in sorted(glob.glob(os.path.join(ROOT, "src/internal/bronstein/*.m"))):
        for m in re.finditer(r"^(\w+)::usage", open(path).read(), re.M):
            names[m.group(1)] = os.path.relpath(path, ROOT)
    for n in TRAGER_NAMES:
        names[n] = "src/internal/trager/Field.m"
    for n in C_NAMES:
        names[n] = "C builtin"

    tests = ""
    for g in TEST_GLOBS:
        for path in glob.glob(os.path.join(ROOT, g)):
            tests += strip_c_comments(open(path).read()) + "\n"

    def covered(n):
        return re.search(r"(?<![\w`$])" + re.escape(n) + r"(?![\w`])", tests) is not None

    missing = [n for n in sorted(names) if n not in EXEMPT and not covered(n)]
    stale = [n for n in sorted(EXEMPT) if n not in names or covered(n)]

    for n in missing:
        print("UNTESTED: %s (%s): add a direct call in a tests/test_bronstein_*.c suite" % (n, names[n]))
    for n in stale:
        print("STALE EXEMPT: %s: %s; delete the entry" %
              (n, "no longer defined" if n not in names else "now tested"))
    print("bronstein coverage: %d routines, %d tested, %d exempt, %d untested, %d stale exemptions"
          % (len(names), len(names) - len(missing) - len(EXEMPT), len(EXEMPT), len(missing), len(stale)))
    return 1 if (missing or stale) else 0


if __name__ == "__main__":
    sys.exit(main())
