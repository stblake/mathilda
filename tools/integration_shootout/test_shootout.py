#!/usr/bin/env python3
"""test_shootout.py -- unit tests of the integration-shootout harness.

  * translate: the FullForm parser and the four emitters on fixed inputs, and
    the refusal of unknown heads;
  * translation round trip: the emitted integrand, evaluated by sympy and by
    Maxima at the corpus points, equals the corpus's own 40-digit value f(p)
    (from Mathilda) -- so a mistranslation can never pass as a system's failure;
  * classify_run / verdicts: truth tables, including the hard verdicts
    (WRONG, FALSE_CERT, DISAGREE) on synthetic records.

Run: python3 tools/integration_shootout/test_shootout.py   (exit 1 on failure;
the round trips are skipped for a system that is not installed).
"""
import os
import re
import subprocess
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import adapters as A  # noqa: E402
import shootout as S  # noqa: E402
import translate as T  # noqa: E402

FAIL = []


def eq(got, want, what):
    if got != want:
        FAIL.append(f"{what}: got {got!r}, want {want!r}")


def test_translate():
    t = T.parse("Plus[Times[x, Power[Sec[x], 2]], Power[E, Times[Rational[1, 2], x]], "
                "Power[Plus[1, x], Rational[-1, 3]], Log[x]]")
    eq(T.emit(t, "sympy"), "((x*(sec(x)**2)) + (E**(Rational(1, 2)*x)) + ((1 + x)**Rational(-1, 3)) + log(x))", "sympy")
    eq(T.emit(t, "maxima"), "((x*(sec(x)^2)) + (%e^((1/2)*x)) + ((1 + x)^(-1/3)) + log(x))", "maxima")
    eq(T.emit(t, "fricas"), "((x*(sec(x)^2)) + (%e^((1/2)*x)) + ((1 + x)^(-1/3)) + log(x))", "fricas")
    eq(T.emit(t, "wolfram"), "((x*(Sec[x]^2)) + (E^((1/2)*x)) + ((1 + x)^(-1/3)) + Log[x])", "wolfram")
    eq(T.emit(T.parse("Times[-2, x]"), "maxima"), "((-2)*x)", "negative integer")
    eq(T.parse("f[]"), ("f",), "empty call")
    for bad in ("BesselJ[0, x]", "Plus[x, y]", "Plus[x"):
        try:
            T.emit(T.parse(bad), "sympy")
            FAIL.append(f"untranslatable accepted: {bad}")
        except (T.Untranslatable, IndexError):
            pass
    eq(T.special_q("sympy", "erf(x) + x"), True, "special sympy")
    eq(T.special_q("sympy", "exp(x) + log(x)"), False, "elementary sympy")
    eq(T.special_q("mathilda", "ExpIntegralEi[x]"), True, "special mathilda")
    eq(T.special_q("maxima", "gamma_incomplete(0,x)"), True, "special maxima")
    eq(T.uneval_q("fricas", "integral(f, x)"), True, "uneval fricas")
    eq(T.uneval_q("sympy", "Integral(1/log(x), x)"), True, "uneval sympy")
    eq(T.uneval_q("maxima", "log(x)"), False, "evaluated maxima")
    # every corpus case translates for every system
    bad = 0
    for c in S.load_corpus():
        for s in T.SYSTEMS:
            try:
                T.emit(c["tree"], s)
            except T.Untranslatable:
                bad += 1
    eq(bad, 0, "untranslatable corpus cases")
    # benchmark 257 too, when it has been generated
    c257 = os.path.join(S.A.REPO, "benchmarks", "257-compact-integrals", "corpus.tsv")
    if os.path.exists(c257):
        bad = sum(1 for c in S.load_corpus(c257) for s in T.SYSTEMS
                  if not _translates(c["tree"], s))
        eq(bad, 0, "untranslatable corpus-257 cases")


def _translates(tree, system):
    try:
        T.emit(tree, system)
        return True
    except T.Untranslatable:
        return False


def corpus_sample(n=40):
    cs = S.load_corpus()
    step = max(1, len(cs) // n)
    return cs[::step][:n]


def fvals(c):
    # the corpus's 40-digit f(p), from column 9 of corpus.tsv
    line = next(ln for ln in open(S.CORPUS) if ln.split("\t")[1] == c["id"])
    vals = line.split("\t")[8].strip("{}").split(",")
    return [Fraction(re.sub(r"`.*", "", v.strip()).replace("*^", "e")) for v in vals]


def close(a, b):
    return abs(a - b) <= Fraction(1, 10 ** 25) * max(1, abs(b))


def test_roundtrip_sympy():
    if not A.available("sympy"):
        print("note: sympy not installed, round trip skipped")
        return
    sample = corpus_sample()
    prog = ["import sympy as sp", "from sympy import *", "x = sp.Symbol('x')"]
    for c in sample:
        for p in c["points_list"]:
            prog.append(f"print({c['id']!r}, sp.N(({T.emit(c['tree'], 'sympy')}).subs(x, sp.Rational({p!r})), 40))")
    out = subprocess.run([A.SYMPY_PY, "-c", "\n".join(prog)], capture_output=True, text=True, timeout=600).stdout
    got = {}
    for ln in out.splitlines():
        cid, v = ln.split(" ", 1)
        got.setdefault(cid, []).append(v.strip())
    for c in sample:
        want = fvals(c)
        vs = got.get(c["id"], [])
        if len(vs) != len(want) or not all(close(Fraction(v), w) for v, w in zip(vs, want)):
            FAIL.append(f"sympy round trip {c['id']}: {vs} vs {[float(w) for w in want]}")


def test_roundtrip_maxima():
    if not A.available("maxima"):
        print("note: maxima not installed, round trip skipped")
        return
    sample = corpus_sample()
    lines = ["display2d:false$ fpprec:45$"]
    for c in sample:
        for p in c["points_list"]:
            lines.append(f'print("@@V {c["id"]}", bfloat(subst(x = {p}, {T.emit(c["tree"], "maxima")})))$')
    out = subprocess.run([A.MAXIMA, "--very-quiet", "--batch-string=" + "\n".join(lines)],
                         capture_output=True, text=True, timeout=600, input="").stdout
    got = {}
    for m in re.finditer(r"^@@V (\S+) (\S+)\s*$", out, re.M):
        got.setdefault(m.group(1), []).append(m.group(2))
    for c in sample:
        want = fvals(c)
        vs = got.get(c["id"], [])
        ok = len(vs) == len(want)
        if ok:
            try:
                ok = all(close(Fraction(v.replace("b", "e")), w) for v, w in zip(vs, want))
            except ValueError:
                ok = False
        if not ok:
            FAIL.append(f"maxima round trip {c['id']}: {vs} vs {[float(w) for w in want]}")


def rec(status, nums=(), cert=False, text=""):
    return {"status": status, "nums": list(nums), "cert": cert, "text": text}


def test_classify():
    eq(S.classify_run("mathilda", rec("ANSWER", ["ok", "ok", "ok"], text="Log[x]")), "ELEM", "elem")
    eq(S.classify_run("mathilda", rec("ANSWER", ["ok", "na", "ok"], text="Erf[x]")), "SPECIAL", "special")
    eq(S.classify_run("sympy", rec("ANSWER", ["ok", "bad", "ok"], text="log(x)")), "WRONG", "wrong")
    eq(S.classify_run("fricas", rec("ANSWER", ["na", "na", "na"])), "UNVERIFIED", "unverified")
    eq(S.classify_run("mathilda", rec("UNEVAL", cert=True)), "CERT", "cert")
    eq(S.classify_run("mathilda", rec("UNEVAL")), "UNEVAL", "uneval")
    eq(S.classify_run("maxima", rec("TIMEOUT")), "TIMEOUT", "timeout")


def test_verdicts():
    elem = {"expect": "elem"}
    unk = {"expect": "unknown"}
    non = {"expect": "nonelem"}
    v = S.verdicts(elem, {"Mathilda:Automatic": "CERT", "FriCAS": "UNEVAL"})
    eq(any(x.startswith("FALSE_CERT") for x in v), True, "false cert")
    v = S.verdicts(unk, {"Mathilda:Bronstein": "CERT", "Mathilda:Automatic": "UNEVAL", "sympy": "ELEM"})
    eq(any(x.startswith("DISAGREE") for x in v), True, "disagree")
    v = S.verdicts(unk, {"Mathilda:Bronstein": "CERT", "sympy": "SPECIAL"})
    eq([x for x in v if x.split()[0] in S.HARD], [], "special answer does not contradict a certificate")
    v = S.verdicts(unk, {"Mathilda:Automatic": "WRONG"})
    eq(any(x.startswith("WRONG") for x in v), True, "wrong")
    v = S.verdicts(unk, {"Mathilda:Automatic": "UNEVAL", "Mathilda:Trager": "ELEM", "FriCAS": "UNEVAL"})
    eq("AUTO-GAP Mathilda:Trager" in v and "ONLY-Mathilda" in v, True, "auto gap + only")
    v = S.verdicts(non, {"Mathilda:Automatic": "ELEM"})
    eq(any(x.startswith("CORPUS?") for x in v), True, "theorem mislabel flagged")
    v = S.verdicts(elem, {"Mathilda:Automatic": "ELEM", "FriCAS": "ELEM", "sympy": "WRONG"})
    eq([x for x in v if x.split()[0] in S.HARD], [], "another system's WRONG is not ours")


def main():
    for t in (test_translate, test_classify, test_verdicts, test_roundtrip_sympy, test_roundtrip_maxima):
        t()
    for f in FAIL:
        print("FAIL", f)
    print("integration shootout harness tests:", "FAILED" if FAIL else "all passed")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
