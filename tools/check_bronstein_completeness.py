#!/usr/bin/env python3
"""check_bronstein_completeness.py -- the completeness ratchet of the Bronstein
completion campaign (BRONSTEIN_PLAN.md, "2026-11-14 -- Completion campaign").

Goal of the campaign: for every integrand in ONE simple radical extension of an
elementary Log/Exp tower, Integrate[f, x, Method -> "Bronstein"] ends in a
verified antiderivative or a certified non-elementarity proof -- never in
"implementation incomplete".  Two lists measure the distance to that goal:

  * every throwUnimplemented["Site:tag", ...] in src/internal/bronstein/*.m and
    src/internal/trager/*.m: an in-scope branch of the algorithm that is not
    written yet (the user sees Integrate::incomplete);
  * every rde_cert_mark("reason") in the C Risch DE (src/calculus/
    integrate_risch_rde.c, risch_field_integrate.c): a step of Bronstein's
    Symbolic Integration I, Ch. 6-7, that cannot certify its "no solution", so
    every decision built on it degrades to Unimplemented.

Each list is diffed against a checked-in BASELINE.  A NEW entry fails (a branch
was added without its implementation, or a new uncertain path appeared); an
entry in BASELINE that no longer exists is reported so the line can be deleted
(the ratchet only moves down).  The campaign is done when both BASELINEs are
empty, the throwUnimplemented class is deleted, and every rde_cert_mark that
remains marks a genuine CAS failure listed in CAS_FAILURE_MARKS.

Usage:  python3 tools/check_bronstein_completeness.py      (exit 1 on a new gap)
"""
import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

M_GLOBS = ["src/internal/bronstein/*.m", "src/internal/trager/*.m"]
C_FILES = ["src/calculus/integrate_risch_rde.c", "src/calculus/risch_field_integrate.c"]

# In-scope branches not implemented yet (site ids of throwUnimplemented).
BASELINE_M = {
    "AlgebraicPart:undecided",       # Risch DE undecided -> Phase 2 (C Ch.6-7)
    "BaseAlgebraic:curve-not-reduced",  # constant-field extension -> Phase 4
    "BaseCurve:curve-not-reduced",   # constant-field extension -> Phase 4
    "Curve:radical-content",         # radical constant field of a component -> Phase 5
    "Curve:radical-of-k",            # a radical below theta -> Phase 5
    "Divisors:multi-ramified",       # d = gcd(e, m) places over a special point -> Phase 3
    "LimitedIntegrate:logrel",       # undecided log relation -> Phase 0/2 structure theorem
    "LogPart:special-genus",         # special residues, genus >= 1 -> Phase 3 (#9)
    "LogPart:special-radical",       # special residue over a radical extension of k -> Phase 3
    "LogPart:weights",               # m >= 3, weight >= 2 coordinates -> Phase 3 (#12)
    "Pipeline:transcendental",       # RischTranscendental undecided -> Phase 2
    "RDE:base-undecided",            # Risch DE over the base undecided -> Phase 2
    "RDEAlgebraic:undecided",        # Risch DE on the curve -> Phase 5
    "Residues:infinity-limited",     # limited integration undecided -> Phase 2/3
    "Residues:nonrational",          # residue polynomial over a larger field -> Phase 4
    "Residues:radical-value",        # special value over a radical extension of k -> Phase 5
    "Residues:special-value",        # value integral undecided -> Phase 2/3
    "Special:rde",                   # special-place Risch DE undecided -> Phase 2
    "Tower:algebraic-constant",      # number-field constants -> Phase 4
    "Tower:complex",                 # Q(i) constants -> Phase 4
    "Tower:compositum-undecided",    # Kummer group not decided -> Phase 3
    "Tower:curve-monomial-tower",    # monomials over a radical-argument kernel -> Phase 5
    "Tower:cyclic-compositum",       # one radical in disguise -> Phase 3
    "Tower:exp-of-log-radical",      # E^(r Log a) as a radical -> Phase 3
    "Tower:exprel",                  # undecided exponential relation -> Phase 3
    "Tower:logrel",                  # undecided logarithm relation -> Phase 3
    "Tower:not-algebraic",           # irrational exponents, parameters -> Phase 4
    "Tower:radical-in-kernel",       # radical inside a kernel above others -> Phase 5
    "Tower:radical-mid-tower",       # radical in the middle of the tower -> Phase 5
    "Tower:radicand-not-rational",   # symbolic parameters -> Phase 4
    "Tower:reducible-log",           # reducible logarithm -> Phase 3 (#13)
    "Tower:trig",                    # TrigToExp front end over Q(i) -> Phase 4
    "Trager",                        # a Trager decline inside Bronstein -> Trager's own sites
    "Trager:algebraic-constant",     # number-field constants -> Phase 4
    "Trager:constants",              # number-field / complex constants -> Phase 4
    "Trager:logpart",                # no log-part construction applies -> Phase 3
    "Trager:not-algebraic",          # parameters / irrational exponents -> Phase 4
    "Trager:radicand-not-rational",  # symbolic parameters -> Phase 4
    "Trager:reducible-radical",      # y^n - g reducible -> Phase 3
    "TorsionModP:frame-root",        # no frame root mod p: try another prime -> Phase 3
}

# rde_cert_mark reasons that are mathematical gaps of the C Risch DE (Phase 2).
BASELINE_C = {
    # an answer of the field integrator outside K (a new logarithm) is not yet
    # reduced to K by the structure theorem -> ParamRischDE / LimitedIntegrate
    "PolyRischDECancel: log-derivative branch integral undecided",
    "RdeBoundDegreePrim: limited integration undecided",
    # Log/Exp monomials that are algebraically dependent (E^(x + x^2) beside E^x,
    # E^(x^2)): the tower builder must rewrite them into an independent tower
    "tower: dependent monomials",
    "RdeBoundDegreeExp: lc(a) = 0",
    "RdeSpecialDenomExp: a(0) = 0",
    "WeakNormalizer: zero",
    "base: identity gate",
    "base: reduced equation not polynomial",
    "tower: antidifferentiation D H = c",
    "tower: identity gate",
    "tower: reduced equation not polynomial",
    "tower: tangent top",
}

# rde_cert_mark reasons that are permanent: a computation the CAS could not do,
# or an input outside the solver's scope, not a missing branch of the algorithm.
CAS_FAILURE_MARKS = {
    # the decision is complete; only a factorisation / linear solve that did not
    # evaluate, or a residue beyond the machine range, leaves it undecided
    "WeakNormalizer: residues not decided",
    "tower: WeakNormalizer residues not decided",
    "RdeSpecialDenomExp: parametric log-derivative undecided",
    "RdeBoundDegreeExp: parametric log-derivative undecided",
    "RdeBoundDegreePrim: log-derivative undecided",
    "RdeBoundDegreePrim: limited integration out of range",
    "PolyRischDECancel: log-derivative test undecided",
    "SPDE: extended gcd unavailable",
    "SPDE(field): extended gcd unavailable",
    "tower: could not be built",
    "tower: coefficients not rational in the tower",
    "tower: kernels left after substitution",
}


def strip_m_comments(text):
    # (* ... *) comments, possibly nested; keep string literals intact
    out, depth, i, in_str = [], 0, 0, False
    while i < len(text):
        c = text[i]
        if depth == 0 and c == '"' and (i == 0 or text[i - 1] != "\\"):
            in_str = not in_str
        if not in_str and text.startswith("(*", i):
            depth += 1; i += 2; continue
        if not in_str and depth and text.startswith("*)", i):
            depth -= 1; i += 2; continue
        if depth == 0:
            out.append(c)
        i += 1
    return "".join(out)


def strip_c_comments(text):
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return re.sub(r"//[^\n]*", " ", text)


def m_sites():
    sites = {}
    for g in M_GLOBS:
        for path in sorted(glob.glob(os.path.join(ROOT, g))):
            src = strip_m_comments(open(path, encoding="utf-8").read())
            for m in re.finditer(r'throwUnimplemented\[\s*"([^"]+)"', src):
                line = src.count("\n", 0, m.start()) + 1
                sites.setdefault(m.group(1), []).append(
                    "%s:%d" % (os.path.relpath(path, ROOT), line))
    return sites


def c_marks():
    marks = {}
    for rel in C_FILES:
        path = os.path.join(ROOT, rel)
        src = strip_c_comments(open(path, encoding="utf-8").read())
        for m in re.finditer(r'rde_cert_mark\(\s*"([^"]+)"\s*\)', src):
            line = src.count("\n", 0, m.start()) + 1
            marks.setdefault(m.group(1), []).append("%s:%d" % (rel, line))
    return marks


def main():
    bad = False
    sites = m_sites()
    marks = c_marks()
    new_m = sorted(set(sites) - BASELINE_M)
    gone_m = sorted(BASELINE_M - set(sites))
    c_known = BASELINE_C | CAS_FAILURE_MARKS
    new_c = sorted(set(marks) - c_known)
    gone_c = sorted(BASELINE_C - set(marks))
    for s in new_m:
        bad = True
        print("NEW unimplemented site %-34s at %s" % (s, ", ".join(sites[s])))
    for s in new_c:
        bad = True
        print("NEW uncertain Risch DE path %-40r at %s" % (s, ", ".join(marks[s])))
    for s in gone_m:
        print("fixed: %s -- delete it from BASELINE_M" % s)
    for s in gone_c:
        print("fixed: %r -- delete it from BASELINE_C" % s)
    print("completeness: %d in-scope unimplemented sites (%d places), %d uncertain "
          "Risch DE paths; baseline %d + %d"
          % (len(sites), sum(len(v) for v in sites.values()),
             len(set(marks) - CAS_FAILURE_MARKS), len(BASELINE_M), len(BASELINE_C)))
    if bad:
        print("FAIL: implement the branch, or (only with a written reason) add it to the baseline")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
