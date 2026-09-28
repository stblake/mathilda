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
    "src/bitwise/bitlength.c:BitLength::argx",
    "src/bitwise/bitlength.c:BitLength::int",
    "src/calculus/integrate.c:Integrate::method",
    "src/calculus/integrate.c:Integrate::nonelem",
    "src/calculus/integrate_line.c:Integrate::idiv",
    "src/calculus/integrate_newton_leibniz.c:Integrate::idiv",
    "src/calculus/limit.c:Limit::method",
    "src/calculus/residue.c:Residue::argm",
    "src/complex.c:Conjugate::argx",
    "src/complex_expand.c:General::argct",
    "src/context.c:Begin::cxt",
    "src/context.c:BeginPackage::cxt",
    "src/context.c:Context::notfound",
    "src/context.c:End::noctx",
    "src/context.c:EndPackage::noctx",
    "src/core.c:%s::rvalue",
    "src/core.c:Negative::argx",
    "src/core.c:NonNegative::argx",
    "src/core.c:NonPositive::argx",
    "src/core.c:Positive::argx",
    "src/core.c:Symbol::symname",
    "src/core.c:Unset::wrsym",
    "src/eval.c:$RecursionLimit::limset",
    "src/eval.c:$RecursionLimit::reclim",
    "src/eval.c:%s::flagset",
    "src/eval.c:%s::nofwd",
    "src/eval.c:%s::wrsym",
    "src/eval.c:Goto::nolabel",
    "src/eval.c:Throw::nocatch",
    "src/expand.c:ExpandAll::argt",
    "src/expand_power.c:PowerExpand::argt",
    "src/facint.c:FactorInteger::nofac",
    "src/fit.c:%s::<dynamic>",
    "src/funcprog.c:MapIndexed::nonopt",
    "src/funcprog.c:MapIndexed::optx",
    "src/int.c:DigitCount::argb",
    "src/int.c:DigitCount::base",
    "src/int.c:DigitCount::digit",
    "src/int.c:DigitCount::int",
    "src/int.c:DigitCount::ovfl",
    "src/int.c:DigitSum::argb",
    "src/int.c:DigitSum::base",
    "src/int.c:DigitSum::int",
    "src/int.c:FromDigits::argb",
    "src/int.c:FromDigits::char",
    "src/int.c:FromDigits::ibase",
    "src/int.c:FromDigits::nlst",
    "src/int.c:IntegerDigits::argb",
    "src/int.c:IntegerDigits::ibase",
    "src/int.c:IntegerDigits::int",
    "src/int.c:IntegerDigits::intnn",
    "src/int.c:IntegerExponent::argt",
    "src/int.c:IntegerExponent::ibase",
    "src/int.c:IntegerExponent::int",
    "src/int.c:IntegerLength::argt",
    "src/int.c:IntegerLength::ibase",
    "src/int.c:IntegerLength::int",
    "src/int.c:IntegerString::argb",
    "src/int.c:IntegerString::basf",
    "src/int.c:IntegerString::int",
    "src/int.c:IntegerString::intnn",
    "src/interp.c:InterpolatingFunction::dmval",
    "src/interp.c:InterpolatingPolynomial::noipf",
    "src/interp.c:InterpolatingPolynomial::poised",
    "src/list/array_pad.c:ArrayPad::mindimsize",
    "src/list/matrixq.c:DiagonalMatrixQ::argt",
    "src/list/matrixq.c:DiagonalMatrixQ::nonopt",
    "src/list/matrixq.c:SquareMatrixQ::argx",
    "src/list/matrixq.c:UpperTriangularMatrixQ::argt",
    "src/list/matrixq.c:UpperTriangularMatrixQ::nonopt",
    "src/list/pad.c:%s::argb",
    "src/list/rescale.c:Rescale::argb",
    "src/loadmodule.c:LoadModule::nofile",
    "src/names.c:Names::regavail",
    "src/names.c:Names::regex",
    "src/nc_accuracy.c:%s::accgl",
    "src/numberform.c:NumberForm::reqsigz",
    "src/numbertheory/divisible.c:Divisible::argm",
    "src/numbertheory/divisible.c:Divisible::argt",
    "src/numbertheory/divisors.c:Divisors::argx",
    "src/numbertheory/divisorsigma.c:DivisorSigma::argrx",
    "src/numbertheory/extendedgcd.c:ExtendedGCD::egcd",
    "src/numbertheory/extendedgcd.c:ExtendedGCD::exact",
    "src/numbertheory/jacobisymbol.c:JacobiSymbol::argrx",
    "src/numbertheory/liouvillelambda.c:LiouvilleLambda::argt",
    "src/numbertheory/moebiusmu.c:MoebiusMu::argx",
    "src/numbertheory/multiplicativeorder.c:MultiplicativeOrder::argt",
    "src/numbertheory/prime.c:Prime::argx",
    "src/numbertheory/prime.c:Prime::intpp",
    "src/numbertheory/prime.c:PrimePi::method",
    "src/numbertheory/primenu.c:PrimeNu::argt",
    "src/numbertheory/primeomega.c:PrimeOmega::argt",
    "src/numbertheory/primitiveroot.c:PrimitiveRoot::argt",
    "src/numbertheory/primitiveroot.c:PrimitiveRoot::intg",
    "src/numbertheory/primitiveroot.c:PrimitiveRootList::argx",
    "src/numeric.c:N::prec",
    "src/options_builtin.c:SetOptions::locked",
    "src/options_builtin.c:SetOptions::optnf",
    "src/partitions.c:IntegerPartitions::argb",
    "src/partitions.c:IntegerPartitions::take",
    "src/partitions.c:IntegerPartitions::undef",
    "src/partitions.c:PartitionsP::argx",
    "src/partitions.c:PartitionsQ::argx",
    "src/precision.c:SetPrecision::prec",
    "src/product/product_infinite.c:Product::div",
    "src/product/product_rational_infinite.c:Product::div",
    "src/purefunc.c:Function::slot1",
    "src/purefunc.c:Function::slota",
    "src/real.c:MantissaExponent::argt",
    "src/real.c:MantissaExponent::ibase",
    "src/real.c:MantissaExponent::realx",
    "src/real.c:RealDigits::argb",
    "src/real.c:RealDigits::ibase",
    "src/real.c:RealDigits::int",
    "src/real.c:RealDigits::intnn",
    "src/real.c:RealDigits::nrep",
    "src/real.c:RealDigits::ovfl",
    "src/real.c:RealDigits::realx",
    "src/real.c:RealExponent::argt",
    "src/real.c:RealExponent::ibase",
    "src/real.c:RealExponent::realx",
    "src/refine.c:Refine::argt",
    "src/refine.c:Refine::argx",
    "src/repl_hooks.c:$PreRead::strret",
    "src/root_numeric.c:Root::<dynamic>",
    "src/rootreduce.c:RootReduce::argx",
    "src/rootreduce.c:RootReduce::mtd",
    "src/solve/reduce.c:Reduce::ivar",
    "src/solve/reduce.c:Reduce::optx",
    "src/solve/reduce_companions.c:FindInstance::optx",
    "src/solve/solve.c:Solve::ivar",
    "src/solve/solve.c:Solve::optx",
    "src/solve/solve_common.c:Solve::svars",
    "src/solve/solvealways.c:SolveAlways::argt",
    "src/solve/solvealways.c:SolveAlways::eqf",
    "src/solve/solvealways.c:SolveAlways::ivar",
    "src/solve/solvelinsys.c:Solve::svars",
    "src/solve/solverad.c:Solve::nongen",
    "src/special_functions/airyai.c:AiryAi::argx",
    "src/special_functions/airyai.c:AiryAiPrime::argx",
    "src/special_functions/airybi.c:AiryBi::argx",
    "src/special_functions/airybi.c:AiryBiPrime::argx",
    "src/special_functions/bernoullib.c:BernoulliB::argt",
    "src/special_functions/bessel.c:BesselI::argrx",
    "src/special_functions/bessel.c:BesselJ::argrx",
    "src/special_functions/bessel.c:BesselK::argrx",
    "src/special_functions/bessel.c:BesselY::argrx",
    "src/special_functions/beta.c:Beta::argb",
    "src/special_functions/coshintegral.c:CoshIntegral::argx",
    "src/special_functions/cosintegral.c:CosIntegral::argx",
    "src/special_functions/erf.c:Erf::argt",
    "src/special_functions/erfc.c:Erfc::argx",
    "src/special_functions/erfi.c:Erfi::argx",
    "src/special_functions/eulere.c:EulerE::argt",
    "src/special_functions/expintegralei.c:ExpIntegralEi::argx",
    "src/special_functions/fresnel.c:%s::argx",
    "src/special_functions/gamma.c:Gamma::argt",
    "src/special_functions/harmonicnumber.c:HarmonicNumber::argt",
    "src/special_functions/hurwitzzeta.c:HurwitzZeta::argrx",
    "src/special_functions/inverf.c:InverseErf::argt",
    "src/special_functions/inverfc.c:InverseErfc::argx",
    "src/special_functions/legendre.c:LegendreP::argb",
    "src/special_functions/legendre.c:LegendreQ::argb",
    "src/special_functions/lerchphi.c:LerchPhi::argrx",
    "src/special_functions/lerchphi.c:LerchPhi::nonopt",
    "src/special_functions/loggamma.c:LogGamma::argx",
    "src/special_functions/logintegral.c:LogIntegral::argx",
    "src/special_functions/pochhammer.c:Pochhammer::argrx",
    "src/special_functions/polygamma.c:PolyGamma::argt",
    "src/special_functions/polylog.c:PolyLog::argt",
    "src/special_functions/productlog.c:ProductLog::argt",
    "src/special_functions/sinc.c:Sinc::argx",
    "src/special_functions/sinhintegral.c:SinhIntegral::argx",
    "src/special_functions/sinintegral.c:SinIntegral::argx",
    "src/special_functions/zeta.c:Zeta::argt",
    "src/strings/regex/regex_common.c:%s::regavail",
    "src/strings/regex/regex_common.c:%s::regex",
    "src/strings/regex/regularexpression.c:RegularExpression::regex",
    "src/strings/stringdrop.c:StringDrop::argrx",
    "src/strings/stringinsert.c:StringInsert::argrx",
    "src/strings/stringreplacepart.c:StringReplacePart::argt",
    "src/strings/stringreplacepart.c:StringReplacePart::ovlp",
    "src/strings/stringreverse.c:StringReverse::argx",
    "src/vectoranal.c:%s::chart",
    "src/vectors.c:UnitVector::nonopt",
    "src/vectors.c:UnitVector::wprec",
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
