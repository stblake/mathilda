#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Refine stress suite — a ratcheting, first-principles regression corpus for
Refine[expr, assum] (and the assumption engine it shares with Simplify /
PossibleZeroQ / Element / Assuming).

Unlike tests/test_refine.c (which asserts a fixed set of PASSing cases in C),
this suite holds ALL 159 cases including the ones Refine still cannot reduce,
and ratchets: it fails if a passing case regresses, if a soundness/robustness
probe ever returns a wrong answer / hangs / crashes, OR if a case on the
known-gap BASELINE starts passing (so the baseline is kept honest and the line
gets deleted). This mirrors the check-* audits in the makefile.

Each case is decided by structural equality (SameQ / ===) inside Mathilda, so
printer-form differences never cause a false failure. Every case runs in its
OWN Mathilda process under a hard timeout, so a hang or crash is contained and
attributed rather than fatal to the run.

Verdicts:
  PASS       result === expected
  UNCHANGED  result === the input   (a missing reduction; never a wrong answer)
  DIVERGENT  neither                (potential wrong answer -> hard failure)
  ERROR      Mathilda exited nonzero with no verdict line
  HANG       process exceeded the per-case timeout
  CRASH      process died by a signal

Usage:
  python3 tests/refine_stress/run_refine_stress.py [--timeout S] [--verbose]
Exit status is nonzero on any regression, hard failure, or stale baseline.
"""
import argparse, collections, json, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
BIN  = os.path.join(REPO, "Mathilda")

sys.path.insert(0, HERE)
from corpus import CASES

# Known gaps: cases that currently return the input unchanged (a missing
# reduction, never a wrong answer). Documented in tasks/refine_stress_report.md.
# Ratchet: if one of these starts PASSing, delete it from this set (the suite
# will tell you to); if a case NOT listed here goes UNCHANGED, that is a
# regression and the suite fails.
BASELINE_UNCHANGED = frozenset({
    "trg08",  # Exp[2 Pi I k] -> 1            (integer-linear trig argument)
    "trg09",  # Sin[(2k+1) Pi/2] -> (-1)^k    (integer-linear trig argument)
    "trg10",  # Cos[2 k Pi] -> 1              (integer-linear trig argument)
    "elt09",  # Element[Sqrt[2], Algebraics]  (algebraic-literal recognition)
    "inq16",  # 7-variable entailment         (CAD bails at > 6 vars)
    "dpp03",  # Abs[(x-1)^2] -> (x-1)^2       (deep-positivity: non-strict sign)
    "dpp04",  # Sqrt[(x-1)^2] -> Abs[x-1]     (deep-positivity: compound real base)
    "adv03",  # Sign of an exact algebraic 0  (needs zero-normalized argument)
})

ACT, PASS, UNCH = "@@ACT@@", "@@PASS@@", "@@UNCH@@"


def build_m(c):
    # The accumulator name must NOT occur in any corpus expression, or it
    # self-references (e.g. `r = Refine[(x^2)^r, ...]`) into the recursion limit.
    a = "RefineStressAcc"
    lines = ["%s = %s;" % (a, c["expr"]),
             'Print["%s", %s];' % (ACT, a),
             'Print["%s", TrueQ[%s === (%s)]];' % (PASS, a, c["expected"])]
    if c.get("baseline") is not None:
        lines.append('Print["%s", TrueQ[%s === (%s)]];' % (UNCH, a, c["baseline"]))
    return "\n".join(lines) + "\n"


def run_case(c, timeout):
    fd, path = tempfile.mkstemp(suffix=".m", dir=HERE)
    with os.fdopen(fd, "w") as f:
        f.write(build_m(c))
    try:
        try:
            p = subprocess.run([BIN, "-file", path], cwd=REPO,
                               capture_output=True, text=True, timeout=timeout)
            rc, out, err = p.returncode, p.stdout, p.stderr
        except subprocess.TimeoutExpired as e:
            return {"id": c["id"], "cat": c["cat"], "verdict": "HANG",
                    "actual": None, "rc": None, "stderr": ""}
    finally:
        os.unlink(path)

    actual = pass_b = unch_b = None
    for ln in out.splitlines():
        if ln.startswith(ACT):    actual = ln[len(ACT):].strip()
        elif ln.startswith(PASS): pass_b = ln[len(PASS):].strip()
        elif ln.startswith(UNCH): unch_b = ln[len(UNCH):].strip()

    if rc < 0:                 verdict = "CRASH"        # killed by a signal
    elif pass_b is None:       verdict = "ERROR"        # no verdict line
    elif pass_b == "True":     verdict = "PASS"
    elif unch_b == "True":     verdict = "UNCHANGED"
    else:                      verdict = "DIVERGENT"
    return {"id": c["id"], "cat": c["cat"], "expr": c["expr"],
            "expected": c["expected"], "note": c["note"], "verdict": verdict,
            "actual": actual, "rc": rc, "stderr": err.strip()[:300]}


def main():
    ap = argparse.ArgumentParser(description="Refine stress suite (ratcheting).")
    ap.add_argument("--timeout", type=float, default=40.0,
                    help="per-case wall-clock timeout in seconds (default 40)")
    ap.add_argument("--verbose", action="store_true",
                    help="print the full category x verdict matrix and every non-PASS case")
    args = ap.parse_args()

    if not os.path.exists(BIN):
        print("ERROR: %s not found. Run `make` first." % BIN, file=sys.stderr)
        return 2

    results = [run_case(c, args.timeout) for c in CASES]
    by_id = {r["id"]: r for r in results}
    counts = collections.Counter(r["verdict"] for r in results)

    if args.verbose:
        verds = ["PASS", "UNCHANGED", "DIVERGENT", "ERROR", "HANG", "CRASH"]
        bycat = collections.OrderedDict()
        for r in results:
            bycat.setdefault(r["cat"], collections.Counter())[r["verdict"]] += 1
        print("%-26s" % "category" + "".join("%-11s" % v for v in verds) + "TOTAL")
        for cat, cnt in bycat.items():
            print("%-26s" % cat[:26] + "".join("%-11d" % cnt.get(v, 0) for v in verds)
                  + "%d" % sum(cnt.values()))
        for r in results:
            if r["verdict"] != "PASS":
                print("  %-10s %-7s %s => %s" % (r["verdict"], r["id"], r["expr"][:56], r["actual"]))

    # Classify against the ratchet baseline.
    hard = [r for r in results if r["verdict"] in ("DIVERGENT", "ERROR", "HANG", "CRASH")]
    regressions = [r for r in results
                   if r["verdict"] == "UNCHANGED" and r["id"] not in BASELINE_UNCHANGED]
    fixed = sorted(i for i in BASELINE_UNCHANGED
                   if i in by_id and by_id[i]["verdict"] == "PASS")
    stale = sorted(i for i in BASELINE_UNCHANGED if i not in by_id)

    print("\nRefine stress: %d cases — %s" % (len(results), dict(counts)))
    print("binary: %s" % BIN)

    ok = True
    if hard:
        ok = False
        print("\nFAIL: %d wrong-answer / hang / crash case(s) — Refine must never do this:" % len(hard))
        for r in hard:
            print("  %-9s %-7s %s\n           expected %s, got %s"
                  % (r["verdict"], r["id"], r["expr"], r["expected"], r["actual"]))
            if r.get("stderr"):
                print("           stderr: %s" % r["stderr"].replace("\n", " | "))
    if regressions:
        ok = False
        print("\nFAIL: %d case(s) newly UNCHANGED (a reduction regressed) — not on the baseline:" % len(regressions))
        for r in regressions:
            print("  %-7s %s\n           expected %s, got %s" % (r["id"], r["expr"], r["expected"], r["actual"]))
    if fixed:
        ok = False
        print("\nFAIL: %d baseline gap(s) now PASS — delete them from BASELINE_UNCHANGED:" % len(fixed))
        print("      " + ", ".join(fixed))
    if stale:
        ok = False
        print("\nFAIL: BASELINE_UNCHANGED names ids not in the corpus: " + ", ".join(stale))

    if ok:
        print("\nAll Refine stress cases at or ahead of baseline "
              "(%d PASS, %d known gaps)." % (counts.get("PASS", 0), len(BASELINE_UNCHANGED)))
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
