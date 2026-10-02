#!/usr/bin/env python3
"""
check_corpus_indvar.py -- does every DSolve corpus record ask the book's question?

THE FAILURE MODE THIS EXISTS FOR. A converted corpus record
(DSolve_test_status/DE_examples_*.m) can PARSE cleanly and still encode a
DIFFERENT equation than the source printed, because the converter inferred the
wrong independent variable. The harness then back-substitutes the solution into
the garbled equation, so the record scores PASS (or UNEVAL) for a question
nobody asked -- a wrong answer that no residual, no verify gate and no gate
baseline can see. Both failure directions are silent.

The recurring shape is a PARAMETER promoted to the independent variable:

    2.2.18-1744   y'' - 2a y' + a^2 y = 0      read as an ODE in `a`
    2.2.33-3245   x'' + k^2 x = 0              read as an ODE in `k`
    2.2.26-2563   y'' + w^2 y = 0              read as an ODE in `w`
    2.2.16-1534   y' = a y^((a-1)/a)           read as dy/da -- and this one sat
                                               in README.md as an unexplained
                                               residue for eight waves

Seventeen records across nine sections were in that state when this script was
written (M63), four of them because nobody re-ran the converter after a fix that
had already repaired them. Hand-grepping found the class six times (M42, M44,
M48, M49, M61, M63) and each new section re-opened the door, so the check is
mechanical now.

TWO RULES, both read off the record itself:

  1. `_missing_x` in the CAS classification is upstream stating that the equation
     does not contain its independent variable (Maple's tag for an autonomous
     ODE). So the indvar must NOT occur free in the equation body. Measured over
     the 495 `_missing_x` records in the committed corpora: 8 violations, all 8
     mis-transcriptions, zero false positives.

  2. The independent variable must be a letter that is a variable by convention
     (INDVAR_OK below) or a Greek variable name. Every other letter -- a, b, c,
     k, m, n, w -- is a parameter in this corpus, and a record carrying one as
     its indvar is the step-3 lone-letter bug in tools/latex_ode_to_mathilda.py.

A record may violate a rule deliberately; it goes in EXEMPT below WITH A REASON.
An unexplained violation is an error. Both rules are assert-empty today.

Usage:  python3 tools/check_corpus_indvar.py      (exit 1 on a finding)
        make check-corpus-indvar
"""

import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CORPORA = os.path.join(ROOT, "DSolve_test_status", "DE_examples_*.m")

# Letters that name a variable in this corpus: the converter's own INDVAR_PREF
# (tools/latex_ode_to_mathilda.py) plus `y`, which is the independent variable of
# a swapped-variable record `x = x(y)` (2.2.1-98/99/100, 2.2.30-2961/2962/2966).
INDVAR_OK = set("x t z s r u v y".split())
# Greek variable names as the converter spells them (its GREEK map's values),
# minus the constant Pi. A Greek letter is a legitimate indvar only for a
# first-order ODE -- `dr/dtheta` (2.2.30-2972/2973/2992) -- which the converter
# gates; here we only accept the spelling.
GREEK_OK = set("nu lam mu alpha beta gam delta om sig theta phi rho kap eps "
               "tau xi zeta eta chi psi".split())

# label -> reason. A deliberate violation, one line each.
EXEMPT = {}

# {"label", eqn, fn, iv, "classif", sympy}  -- fn is a symbol or a {list} (system)
REC = re.compile(
    r'^\s*\{"(?P<label>[^"]+)",\s*(?P<eqn>.*),\s*'
    r'(?P<fn>\{[^{}]*\}|[A-Za-z][A-Za-z0-9]*),\s*'
    r'(?P<iv>[A-Za-z][A-Za-z0-9]*),\s*"(?P<cl>.*)",\s*(?:True|False)\}'
)


def iv_free_in_body(eqn, iv):
    """True when `iv` occurs in the equation OUTSIDE every f[iv] application.

    The independent variable is always present as the argument of the dependent
    function, so that spelling has to be blanked before the scan or every record
    looks like a violation.
    """
    stripped = re.sub(r"[A-Za-z][A-Za-z0-9]*'*\[\s*" + re.escape(iv) + r"\s*\]",
                      "@", eqn)
    return re.search(r"(?<![A-Za-z0-9])" + re.escape(iv) + r"(?![A-Za-z0-9])",
                     stripped) is not None


def main():
    # Optional explicit paths so the check can be pointed at a candidate corpus
    # before it is committed, or at a pristine copy to show it catches a known
    # bad state. No argument audits every committed corpus.
    paths = sys.argv[1:] or sorted(glob.glob(CORPORA))
    problems = []
    n_rec = n_missing_x = 0
    unparsed = []
    for path in paths:
        rel = os.path.relpath(path, ROOT)
        for lineno, line in enumerate(open(path, encoding="utf-8"), 1):
            if not line.lstrip().startswith('{"'):
                continue
            m = REC.match(line.rstrip().rstrip(","))
            if not m:
                unparsed.append("%s:%d" % (rel, lineno))
                continue
            n_rec += 1
            label, iv, cl, eqn = (m.group("label"), m.group("iv"),
                                  m.group("cl"), m.group("eqn"))
            if label in EXEMPT:
                continue
            if "_missing_x" in cl:
                n_missing_x += 1
                if iv_free_in_body(eqn, iv):
                    problems.append(
                        "%s  %s: classification says _missing_x, but the "
                        "independent variable `%s` occurs in the equation -- so "
                        "`%s` is a PARAMETER and this record is a different ODE "
                        "than the source printed.\n      %s"
                        % (rel, label, iv, iv, eqn.strip()[:110]))
            if iv not in INDVAR_OK and iv not in GREEK_OK:
                problems.append(
                    "%s  %s: independent variable `%s` is a parameter letter, "
                    "not a variable -- the lone-letter adoption in "
                    "tools/latex_ode_to_mathilda.py misread it.\n      %s"
                    % (rel, label, iv, eqn.strip()[:110]))

    print("read %d records (%d classified _missing_x) from %d file(s)"
          % (n_rec, n_missing_x, len(paths)))
    if unparsed:
        print("  note: %d line(s) did not match the record shape: %s"
              % (len(unparsed), ", ".join(unparsed[:5])))
    for label, why in sorted(EXEMPT.items()):
        print("  exempt: %s -- %s" % (label, why))

    if problems:
        print("\ncorpus independent-variable audit FAILED:\n", file=sys.stderr)
        for p in problems:
            print("  " + p + "\n", file=sys.stderr)
        print("%d problem(s). Re-run tools/latex_ode_to_mathilda.py on the "
              "section's Ch2.S1.SSN.htm page and patch the record, or add an "
              "EXEMPT entry with a reason. See DSolve_test_status/README.md."
              % len(problems), file=sys.stderr)
        return 1

    print("OK: every record's independent variable is a variable, and no "
          "_missing_x record names it in the equation.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
