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
C_FILES = ["src/calculus/integrate_risch_rde.c", "src/calculus/risch_field_integrate.c",
           # the decision's own UNKNOWN exits (marked since 2026-11-19, so a
           # silent decline of Risch`ElementaryIntegralQ is counted here too)
           "src/calculus/integrate_risch_transcendental.c"]

# In-scope branches not implemented yet (site ids of throwUnimplemented).
BASELINE_M = {
    "AlgebraicPart:undecided",       # Risch DE undecided -> Phase 2 (C Ch.6-7)
    "BaseAlgebraic:curve-not-reduced",  # constant-field extension -> Phase 4
    "BaseCurve:curve-not-reduced",   # constant-field extension -> Phase 4
    "Curve:radical-content",         # radical constant field of a component -> Phase 5
    "Curve:radical-of-k",            # a radical below theta -> Phase 5
    "Divisors:multi-ramified",       # d = gcd(e, m) places over a special point -> Phase 3
    "LimitedIntegrate:logrel",       # undecided log relation -> Phase 0/2 structure theorem
    "LogPart:special-radical",       # special residue over a radical extension of k -> Phase 3
    "LogPart:weights",               # m >= 3, weight >= 2 coordinates -> Phase 3 (#12)
    "Pipeline:transcendental",       # RischTranscendental undecided -> Phase 2
    "RDE:base-undecided",            # Risch DE over the base undecided -> Phase 2
    "RDEAlgebraic:undecided",        # Risch DE on the curve -> Phase 5
    "Residues:gaussian-kappa",       # special residues with kappa^(1/m), kappa in Q(i) \ Q -> Phase 4.1
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
    "Tower:exprel",                  # undecided exponential relation -> Phase 3
    "Tower:logrel",                  # undecided logarithm relation -> Phase 3
    "Tower:not-algebraic",           # irrational exponents (and radicals of parameters, Sqrt[a]) -> Phase 4
    "Tower:radical-in-kernel",       # radical inside a kernel above others -> Phase 5
    "Tower:radical-mid-tower",       # radical in the middle of the tower -> Phase 5
    "Tower:radicand-not-rational",   # a radicand with a radical of a parameter (Sqrt[Sqrt[a] + Log[x]]): Q(a)-bar constants -> Phase 4.4 later
    "Trager:parameters-with-constants",  # parameters with Q(i) / radical constants: Phase 4.4 later
    "Bronstein:parameters-with-constants",  # the same constant field Q(i, constants)(params) in Bronstein -> Phase 4.4 later
    "Pipeline:parameter-field-value",  # a log term over Q(params)[t]/(mu), mu of degree >= 3 not binomial, outside the per-orbit construction (no root in radicals to present it) -> Phase 4.4 later
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
    "TorsionModP:complex",           # Q(i) data mod p: map I to a square root of -1 mod p -> Phase 4.1
    "TorsionModP:algebraic-constant",  # radical constants mod p (a degree-1 prime of the constant field) -> Phase 4.3b
    "Trager:constant-outside-field",   # a radical constant formed outside the residue field -> Phase 4.3b
    "Trager:residue-field-QI",       # residue field over Q(i) past ToNumberField (tower method) -> Phase 4.1
    "Trager:residue-field-NF",       # residue field over Q(radical constants) past ToNumberField: the tower method ignores the constants (review of 4.3b) -> G-perf-4
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
    "base: identity gate",
    "base: reduced equation not polynomial",
    "tower: antidifferentiation D H = c",
    "tower: hyperbolic-tangent top",
    # a theta^(>=1) coefficient whose integral carries an ArcTan over a tower
    # with a complex logarithm, and the derivative fold (ArcTan' = sum c_j t_j')
    # did not absorb it: the dilogarithm obstruction (Thm 5.6.1(ii)) is proved
    # only for real towers (batch E: it certified D[Log[(u - I)/(u + I)]^2] not
    # elementary); the general case needs the structure theorem over Q(i)
    "ArcTan over a complex tower: not folded",
    # Phase 2b F13 (review 2026-11-16): limited integration used to decline
    # SILENTLY here, invisible to this gate.  Now each undecided step is marked:
    # - a new logarithm outside rt_newlog_certificate's proof (F5): a log whose
    #   argument holds a tower variable (Log[1 + Log[u]]), an exponential tower,
    #   a non-constant log coefficient (Log[u] Log[1 + u]), parameters in a
    #   non-rational position;
    "LimitedIntegrate: new logarithm undecided",
    # - the lower-field integral itself undecided (no certificate raised);
    "LimitedIntegrate: lower integral undecided",
    # - an integral of the coefficient outside K + C t (degree >= 2 in t_L).
    "LimitedIntegrate: coefficient integral outside K + C t",
    # Phase 2b F8/F13 (v0.6599): the field Risch DE over a lower field whose top
    # is a hyperbolic tangent (sigma = -1, special part t^2 - 1 splits) is a
    # scope decline; the hypertangent top runs rde_tower since v0.6601.
    "RDE: lower field outside the recursive RDE (hyperbolic tangent)",
    # Phase 2b F9 (v0.6599): over algebraic constants the base RDE's polynomial
    # solve / SPDE zero tests are only as good as the normaliser; a failure there
    # is not a proof (a released false certificate: q' + (2 + Sqrt[2]) q = p).
    # Step 7 (v0.6606): the base RDE over numeric algebraic constants runs in an
    # exact mode (RootReduced coefficients, exact division/gcd, exact integer
    # roots).  Outside it -- a parameter beside the constants (Q(a)(Sqrt[2]), the
    # parameters-with-constants scope), or a compositum above degree 8 (RootReduce
    # cost; an AlgebraicNumber representation would lift the cap) -- a "no" is
    # not a theorem.
    "base: algebraic constants outside the exact mode",
    # ... and the tower stack (k(x)[tau] arithmetic) has no exact mode over Q(alpha)
    "tower: algebraic constants (no exact tower arithmetic)",
    # Phase 2b F8 (v0.6601): the hypertangent top.  a(i) = 0 cannot occur after
    # the normal denominator (a is special-free) and a residual not divisible by
    # t^2 + 1 is a normaliser failure; both are defensive.  CancelTan's q(i) must
    # split into real parts in k, and its k(i) Risch DE is unique only when the
    # homogeneous equation has no solution there (Bronstein 6.6: the full
    # CoupledDECancelTan recursion would decide that case).
    "PolyRischDECancelTan: q(i) not split into real parts",
    "PolyRischDECancelTan: homogeneous solution over k(i) not excluded",
    # a tower of several tangents needs the tangent structure theorem to prove
    # the monomials independent (Tan[2 x] = 2 Tan[x]/(1 - Tan[x]^2))
    "tower: more than one tangent",
    # over k(i) t^2 + 1 splits (t - i, t + i): the hypertangent boxes assume it
    # irreducible; the split case is the exponential-like two-factor special part
    "tower: hypertangent top over a field containing I",
    # step 3c: a Liouvillian primitive level whose Int eta is elementary below it
    # (Erfc beside Erf, li[u] beside Ei[Log[u]]) is dependent, and one whose Int eta
    # is undecided is unproved; eliminating a dependent primitive from the tower
    # (the step-8 analogue for logarithms) is the follow-up
    "tower: Liouvillian primitive level not proved a monomial",
    # step 3b: LimitedIntegrate at the bottom coefficient (new logarithms allowed)
    # is decided when eta = g tau^k with tau the exponential directly below theta;
    # li (eta = u'/Log u, a normal pole) and an exponential with other levels
    # between it and theta (Log[x] Erf[1 + x]) are the follow-up
    "LimitedIntegrate(prim): bottom coefficient undecided (li, or an exponential not directly below)",
    # step 3a: the parametric Risch DE does not convert the tangent boxes
    # (NoCancel2/3 and CancelTan's re/im split, where rc_conj takes every symbol,
    # the unknowns too, as real); step 4/5 splits the unknowns into re/im parts
    "parametric: hypertangent top",
    # The decision's UNKNOWN exits (rt_decide_field / rt_decide_tower), silent
    # until 2026-11-19 (32 of the hunts' undecided cases carried no mark at all):
    # Together does not combine a rational function over algebraic constants
    # with symbolic numerators (b/(Sqrt[2] + 2 u) + a/(1 + Sqrt[2] u) stays a
    # sum), so the gate cannot read F as rational in the tower; the fix is the
    # exact tower arithmetic over Q(alpha) (step-7 follow-up)
    "decide: integrand not rational in the tower",
    # kernel shapes the tower builder does not form
    "decide: tower could not be built",
    # no antiderivative and no proof; when it is the only mark, a decline inside
    # the field integrator without its own mark (to be traced case by case)
    "decide: field integrator gave no verdict",
    # the builder's dependence flag (composite exponential classes it does not
    # rewrite); the decision then gives no certificate, as for the RDE
    "decide: dependent monomials",
    # step 3c: Int eta elementary or undecided below a primitive level
    "decide: Liouvillian level not proved a monomial",
    # a constant field beyond Q(i)(params, Pi): Log[2], E, EulerGamma, ... with a
    # tower log that can interact with it (T3, T12)
    "decide: transcendental constant",
    # a hypertangent "no" is a theorem only over k = C(x), real, Tan: over a taller
    # tower the driver solves its coupled system over C(x) only; Tanh and k(i) are
    # the split special part (steps 4, 5)
    "hypertangent top: IntegrateHypertangent's no over a tower, Tanh, or k(i)",
}

# rde_cert_mark reasons that are permanent: a computation the CAS could not do,
# or an input outside the solver's scope, not a missing branch of the algorithm.
CAS_FAILURE_MARKS = {
    # Step 0.3 of the C completeness campaign (v0.6603): defensive marks a hunt
    # over every constant field never reached; each is a state the mathematics rules out.
    # deg a = deg b >= 0 makes lc(a) nonzero; a zero reading is an uncanonical coefficient
    "RdeBoundDegreeExp: lc(a) = 0",
    # a = dn h is built special-free (rc_normal_part strips tau^v), so tau does not divide a
    "RdeSpecialDenomExp: a(0) = 0",
    # q starts at 1 and multiplies gcds of a nonzero d1, never 0
    "WeakNormalizer: zero",
    # a returned y failing the exact identity is a reduction or normaliser slip, not a branch
    "tower: identity gate",
    # a, b, c are polynomials by construction after the normal/special denominators over a real k
    "tower: reduced equation not polynomial",
    # a is special-free and t^2 + 1 is irreducible over a real k (the k(i) case is marked separately)
    "RdeSpecialDenomTan: a(i) = 0",
    # r(i) = 0 by construction and r is real, so t^2 + 1 divides r
    "PolyRischDECancelTan: residual not divisible by tau^2 + 1",
    # the decision is complete; only a factorisation / linear solve that did not
    # evaluate, or a residue beyond the machine range, leaves it undecided
    "WeakNormalizer: residues not decided",
    "tower: WeakNormalizer residues not decided",
    "RdeSpecialDenomExp: parametric log-derivative undecided",
    "RdeBoundDegreeExp: parametric log-derivative undecided",
    "RdeBoundDegreePrim: log-derivative undecided",
    "RdeBoundDegreePrim: limited integration out of range",
    "PolyRischDECancel: log-derivative test undecided",
    # the same decision on the parametric pass (step 3a)
    "parametric: PolyRischDECancel log-derivative test undecided",
    "RdeSpecialDenomTan: parametric log-derivative undecided",
    "RdeBoundDegreeTan: bound out of range",
    "PolyRischDENoCancel3: resonance out of range",
    "SPDE: extended gcd unavailable",
    "SPDE(field): extended gcd unavailable",
    "tower: could not be built",
    "tower: coefficients not rational in the tower",
    "tower: kernels left after substitution",
    # the canonical split of a well-formed tower element does not fail
    "LimitedIntegrate(prim): Laurent coefficient not extracted",
    # outside the transcendental decision's scope: an algebraic function of x is
    # the algebraic integrator's (Trager / Bronstein package), not this one's
    "decide: algebraic function of x",
    # Risch`RischDECertified / ParametricRischDE input gate: f, g outside C(x) with
    # no tower kernel (an input outside the solver's scope)
    "base: coefficients not rational in x",
    # a rewrite, zero test or polynomial primitive that did not evaluate
    "decide: kernels left after substitution",
    "decide: antiderivative not verified",
    "decide: a polynomial primitive did not evaluate",
    # the hypertangent driver is total on a well-formed k(t) element; defensive
    "hypertangent top: the driver did not evaluate",
    "hypertangent top: the driver's remainder is not in k",
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
