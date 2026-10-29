#!/usr/bin/env python3
"""
trager_fricas_stress.py -- nightly differential stress test of Mathilda's
Trager integrator (Integrate[f, x, Method -> "Trager" | "TragerKauers"])
against the Risch--Trager--Bronstein integrator of FriCAS.

Every integrand of tests/data/trager_stress_corpus.m is an algebraic function
over ONE simple radical r(x)^(1/n) with rational coefficients.  Each one is run

  * through Mathilda's raw workers Integrate`Trager[f, x] and
    Integrate`TragerKauers[f, x] (one fresh Mathilda process per case and
    method, the worker under TimeConstrained, plus a process timeout, so a
    crash or hang loses only that one run).  The worker returns the
    antiderivative, {"not elementary", reason} (a proof) or {"failed", reason}
    (a decline).  An antiderivative F is verified inside Mathilda by
    PossibleZeroQ[D[F, x] - f] AND by a 40-digit numeric check of
    D[F, x] - f at 3 rational points;
  * through FriCAS (`fricas -nosman`, one process per case, per-case
    timeout).  Its answer r is verified inside FriCAS by
    complexNumeric(eval(D(r, x) - f, x, p), 40) at the same 3 points.  A
    result containing FriCAS's own unevaluated integral(...) counts as
    "non-elementary" (FriCAS's algebraic integrator is complete for this
    class, so an unevaluated integral is its proof), never as an answer.

The points are chosen by Mathilda (phase 0) from a fixed candidate list so that
f is finite and every radicand is a positive real: derivatives, not
antiderivative values, are compared, and both systems take the principal real
root there, so branch cuts do not enter.

Classification (per Mathilda method):
  AGREE-ELEM          both give a verified antiderivative
  AGREE-NONELEM       Mathilda proves non-elementary, FriCAS stays unevaluated
  MATHILDA-ONLY       verified antiderivative in Mathilda only (when FriCAS
                      stayed unevaluated this is a FriCAS bug: the verified
                      antiderivative is itself the proof of elementarity)
  FRICAS-ONLY         FriCAS gives a verified antiderivative, Mathilda declined
                      ("failed" / unevaluated) or timed out          [ratchet]
  DISAGREE            Mathilda proves non-elementary while FriCAS gives a
                      verified antiderivative, or while the corpus' theorem-
                      backed `expect` field says "elem"            [HARD fail]
  WRONG               a Mathilda antiderivative fails its derivative check
                                                                   [HARD fail]
  NONELEM-UNCONFIRMED Mathilda proves non-elementary, FriCAS errored/timed out
                      or gave a wrong answer
  NEITHER             Mathilda declined, FriCAS has no verified antiderivative
  TIMEOUT / ERROR     the Mathilda run timed out (FriCAS not elementary) / the
                      Mathilda process crashed or printed no result

A FriCAS antiderivative that fails ITS derivative check is reported in the
"FriCAS wrong" list and treated as "no answer" for the comparison (it is not a
Mathilda failure).

Usage (from the repository root):
    python3 tools/trager_fricas_stress.py                 # full run, reports
    python3 tools/trager_fricas_stress.py --check         # nightly gate
    python3 tools/trager_fricas_stress.py --filter goursat --limit 10 --jobs 6
    python3 tools/trager_fricas_stress.py --timeout 30
    python3 tools/trager_fricas_stress.py --emit-random   # regenerate the
        seeded "random" and "deriv" corpus sections (printed to stdout)

--check exits 1 on any WRONG or DISAGREE, or when a method's FRICAS-ONLY count
exceeds BASELINE[method] (lower BASELINE when gaps close -- it ratchets).
Without fricas on PATH (or at /opt/local/bin/fricas) it prints a note and exits
0.  Python 3 standard library only.
"""
import argparse
import concurrent.futures as cf
import json
import os
import random
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import threading
import time
from fractions import Fraction

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CORPUS = os.path.join(REPO, "tests", "data", "trager_stress_corpus.m")
REPORT_MD = os.path.join(REPO, "tests", "data", "trager_stress_report.md")
REPORT_JSON = os.path.join(REPO, "tests", "data", "trager_stress_report.json")

METHODS = ("Trager", "TragerKauers")

# Ratchet: the FRICAS-ONLY count observed on the full corpus.  Lower it when a
# gap closes; --check fails when the count rises above it.
BASELINE = {"Trager": 0, "TragerKauers": 58}

HARD = ("WRONG", "DISAGREE")
CLASSES = ("AGREE-ELEM", "AGREE-NONELEM", "MATHILDA-ONLY", "FRICAS-ONLY",
           "DISAGREE", "WRONG", "NONELEM-UNCONFIRMED", "NEITHER", "TIMEOUT",
           "ERROR")

# Candidate evaluation points (phase 0 keeps the first three good ones).
CANDIDATES = ["3/7", "5/4", "7/3", "13/5", "2/9", "11/3", "17/4", "6/5",
              "29/7", "1/11", "23/3", "31/9", "-1/3", "-3/7", "-5/4", "-7/3",
              "-13/5", "-11/3", "-17/4", "-29/7"]

TOL_EXP = 20            # |D[F] - f| <= 10^-20 * Max[1, |f|]
VERIFY_BUDGET = 60      # seconds of verification on top of --timeout


def find_fricas():
    p = os.environ.get("FRICAS_BIN") or shutil.which("fricas")
    if not p and os.path.exists("/opt/local/bin/fricas"):
        p = "/opt/local/bin/fricas"
    return p


def find_mathilda():
    return os.environ.get("MATHILDA_BIN") or os.path.join(REPO, "Mathilda")


# --------------------------------------------------------------------------
# Corpus
# --------------------------------------------------------------------------

LINE_RE = re.compile(r'^\s*\{\s*"([^"]+)"\s*,\s*"([^"]+)"\s*,\s*"([^"]+)"\s*,'
                     r'\s*"([^"]*)"\s*\}\s*,?\s*$')
TOKEN_RE = re.compile(r"\s*(Sqrt\[|\d+|x|[-+*/^()\]])")


def check_syntax(f):
    """The corpus syntax: integers, x, + - * / ^ ( ), Sqrt[ ].  Products are
    always written with '*' (FriCAS reads juxtaposition as application)."""
    pos, prev = 0, None
    depth = []
    for m in TOKEN_RE.finditer(f):
        if m.start() != pos:
            raise ValueError(f"bad character at {pos}: {f!r}")
        tok = m.group(1)
        if tok in ("x", "(", "Sqrt[") or tok.isdigit():
            if prev is not None and (prev in ("x", ")", "]") or prev.isdigit()):
                raise ValueError(f"implicit product before {tok!r}: {f!r}")
        if tok in ("(", "Sqrt["):
            depth.append(tok)
        elif tok in (")", "]"):
            if not depth or (depth.pop() == "(") != (tok == ")"):
                raise ValueError(f"unbalanced {tok!r}: {f!r}")
        prev = tok
        pos = m.end()
    if pos != len(f.rstrip()) or depth:
        raise ValueError(f"unparsable or unbalanced: {f!r}")


def to_fricas(f):
    check_syntax(f)
    return f.replace("Sqrt[", "sqrt(").replace("]", ")").replace(" ", "")


def load_corpus(path):
    cases, seen = [], set()
    with open(path) as fh:
        text = fh.read()
    # blank out (* ... *) comments, keeping line numbers
    text = re.sub(r"\(\*.*?\*\)", lambda m: "\n" * m.group(0).count("\n"),
                  text, flags=re.S)
    if True:
        for ln, line in enumerate(text.split("\n"), 1):
            m = LINE_RE.match(line)
            if not m:
                continue
            cid, fam, f, expect = m.groups()
            if cid in seen:
                raise SystemExit(f"{path}:{ln}: duplicate id {cid}")
            seen.add(cid)
            check_syntax(f)
            cases.append({"id": cid, "family": fam, "f": f, "expect": expect})
    return cases


# --------------------------------------------------------------------------
# Seeded generators (the "deriv-*" and "random" corpus sections)
# --------------------------------------------------------------------------
# Polynomials are lists of ints, lowest degree first.

def ptrim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p or [0]


def padd(a, b):
    n = max(len(a), len(b))
    return ptrim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                  for i in range(n)])


def pscale(a, c):
    return ptrim([c * v for v in a])


def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            r[i + j] += u * v
    return ptrim(r)


def pder(a):
    return ptrim([i * a[i] for i in range(1, len(a))]) if len(a) > 1 else [0]


def pdeg(a):
    a = ptrim(a)
    return -1 if a == [0] else len(a) - 1


def pgcd_deg(a, b):
    """degree of gcd(a, b) over Q (Euclid on Fractions)."""
    a = [Fraction(v) for v in ptrim(a)]
    b = [Fraction(v) for v in ptrim(b)]
    while pdeg(b) >= 0 and b != [0]:
        while pdeg(a) >= pdeg(b) and a != [0]:
            c = a[-1] / b[-1]
            s = pdeg(a) - pdeg(b)
            a = ptrim([a[i] - (c * b[i - s] if i >= s else 0)
                       for i in range(len(a))])
            if a == [0]:
                break
        a, b = b, a
    return pdeg(a)


def squarefree(p):
    return pdeg(p) >= 1 and pgcd_deg(p, pder(p)) == 0


def pstr(p):
    """Polynomial in corpus syntax, highest degree first, explicit '*'."""
    p = ptrim(p)
    terms = []
    for i in range(len(p) - 1, -1, -1):
        c = p[i]
        if c == 0:
            continue
        mono = "" if i == 0 else ("x" if i == 1 else f"x^{i}")
        a = abs(c)
        if mono == "":
            body = str(a)
        elif a == 1:
            body = mono
        else:
            body = f"{a}*{mono}"
        terms.append(("-" if c < 0 else "+", body))
    if not terms:
        return "0"
    s = ("-" if terms[0][0] == "-" else "") + terms[0][1]
    for sg, body in terms[1:]:
        s += f" {sg} {body}"
    return s


def par(p):
    s = pstr(p)
    return s if re.fullmatch(r"x|\d+", s) else f"({s})"


def rand_poly(rng, deg, lo=-3, hi=3, monic=False):
    while True:
        p = [rng.randint(lo, hi) for _ in range(deg + 1)]
        if monic:
            p[-1] = 1
        if p[-1] != 0:
            return p


def rand_radicand(rng, deg):
    while True:
        r = rand_poly(rng, deg, monic=rng.random() < 0.6)
        if squarefree(r):
            return r


def gen_seeded(seed=20261029):
    """The two seeded sections: 'deriv-alg' / 'deriv-log' are derivatives of
    known elementary functions (expect elem); 'random' is R(x) r(x)^(k/n) with
    small integer coefficients (no expectation)."""
    rng = random.Random(seed)
    out = []
    # d/dx (A r^(k/n)) = (n r A' + k A r') / (n r) * r^(k/n)
    for i in range(20):
        n = rng.choice([2, 2, 3, 3, 4, 5, 6])
        k = rng.choice([kk for kk in range(-n + 1, n) if kk != 0 and
                        Fraction(kk, n).denominator == n])
        r = rand_radicand(rng, rng.randint(2, 5))
        A = rand_poly(rng, rng.randint(0, 2))
        num = padd(pscale(pmul(r, pder(A)), n), pscale(pmul(A, pder(r)), k))
        if num == [0]:
            continue
        e = Fraction(k, n)
        out.append((f"deriv-alg-{i + 1:02d}", "deriv-alg",
                    f"{par(num)}/({n}*{par(r)})*{par(r)}^({e.numerator}/{e.denominator})",
                    "elem"))
    # d/dx (1/2) log((P + Q y)/(P - Q y)), y^2 = r:
    #   y * ((2 r Q' + Q r') P - 2 r P' Q) / (2 r (P^2 - Q^2 r))
    i = 0
    while i < 25:
        r = rand_radicand(rng, rng.randint(3, 6))
        P = rand_poly(rng, rng.randint(1, 3))
        Q = rand_poly(rng, rng.randint(0, 1))
        num = padd(pmul(padd(pscale(pmul(r, pder(Q)), 2), pmul(Q, pder(r))), P),
                   pscale(pmul(pmul(r, pder(P)), Q), -2))
        den = padd(pmul(P, P), pscale(pmul(pmul(Q, Q), r), -1))
        if num == [0] or pdeg(den) < 1:
            continue
        i += 1
        out.append((f"deriv-log-{i:02d}", "deriv-log",
                    f"{par(num)}/(2*{par(r)}*{par(den)})*Sqrt[{pstr(r)}]", "elem"))
    # random R(x) r(x)^(k/n)
    for i in range(60):
        n = rng.choice([2, 2, 2, 3, 3, 3, 4, 4, 5, 6])
        k = rng.choice([kk for kk in range(-n + 1, n) if kk != 0 and
                        Fraction(kk, n).denominator == n])
        r = rand_radicand(rng, rng.randint(1, 4))
        p = rand_poly(rng, rng.randint(0, 2))
        q = rand_poly(rng, rng.randint(0, 2))
        e = Fraction(k, n)
        rad = (f"Sqrt[{pstr(r)}]" if e == Fraction(1, 2)
               else f"{par(r)}^({e.numerator}/{e.denominator})")
        f = f"{par(p)}*{rad}" if pdeg(q) <= 0 and q == [1] else \
            f"{par(p)}/{par(q)}*{rad}"
        out.append((f"random-{i + 1:03d}", "random", f, ""))
    return out


def corpus_line(cid, fam, f, expect):
    return f'  {{"{cid}", "{fam}", "{f}", "{expect}"}}'


# --------------------------------------------------------------------------
# Mathilda
# --------------------------------------------------------------------------

def run_proc(cmd, timeout, stdin=None):
    t0 = time.time()
    try:
        pr = subprocess.run(cmd, input=stdin, capture_output=True, text=True,
                            timeout=timeout, cwd=REPO)
        return pr.stdout, pr.stderr, pr.returncode, time.time() - t0, False
    except subprocess.TimeoutExpired as e:
        def dec(s):
            return s.decode(errors="replace") if isinstance(s, bytes) else (s or "")
        return dec(e.stdout), dec(e.stderr), None, time.time() - t0, True


def mathilda_script(path, body):
    with open(path, "w") as fh:
        fh.write(body)


POINTS_TMPL = r"""
Stress`OkQ[f_, bases_, p_] := Module[{fv, bv},
  fv = Quiet[N[f /. x -> p, 30]]; bv = Quiet[N[bases /. x -> p, 30]];
  NumberQ[fv] && Im[fv] == 0 && Abs[fv] < 10^8 &&
    AllTrue[bv, (NumberQ[#] && Im[#] == 0 && # > 1/100) &]];
Stress`Pts[f_] := Module[{bases, good},
  bases = Cases[{f}, Power[b_, e_Rational] :> b, {0, Infinity}];
  good = Select[CANDS, Stress`OkQ[f, bases, #] &];
  If[Length[good] < 3,
    good = Join[good, Select[CANDS, (!MemberQ[good, #] &&
      NumberQ[Quiet[N[f /. x -> #, 30]]]) &]]];
  Take[good, UpTo[3]]];
"""


def mathilda_points(mbin, cases, tmpdir):
    """Phase 0: three evaluation points per case, all cases in one process
    (with a per-case fallback if that process fails)."""
    cands = "{" + ", ".join(CANDIDATES) + "}"
    pre = POINTS_TMPL.replace("CANDS", cands)

    def batch(cs):
        body = pre + "".join(
            f'Print["@@PTS {c["id"]} ", InputForm[TimeConstrained['
            f'Stress`Pts[{c["f"]}], 20, {{}}]]];\n' for c in cs)
        path = os.path.join(tmpdir, f"pts_{cs[0]['id']}.m")
        mathilda_script(path, body)
        out, _, rc, _, to = run_proc([mbin, "-file", path], 30 + 20 * len(cs))
        res = {}
        for m in re.finditer(r"^@@PTS (\S+) \{(.*)\}\s*$", out, re.M):
            res[m.group(1)] = [s.strip() for s in m.group(2).split(",") if s.strip()]
        return res

    res = {}
    for i in range(0, len(cases), 50):
        res.update(batch(cases[i:i + 50]))
    for c in cases:
        if c["id"] not in res:
            res.update(batch([c]))
        c["points"] = res.get(c["id"]) or CANDIDATES[:3]


RUN_TMPL = r"""
f = (INTEGRAND);
pts = POINTS;
{t, r} = AbsoluteTiming[TimeConstrained[Quiet[WORKER[f, x]], TIMEOUT, $Aborted]];
Print["@@TIME ", InputForm[N[t]]];
Print["@@RES ", InputForm[r]];
status = Which[
  r === $Aborted, "timeout",
  ListQ[r] && Length[r] == 2 && StringQ[r[[1]]], r[[1]],
  !FreeQ[r, Integrate | WORKER], "unevaluated",
  True, "elem"];
Print["@@STATUS ", status];
If[status === "elem",
  d = Quiet[D[r, x] - f];
  Print["@@PZQ ", TimeConstrained[Quiet[PossibleZeroQ[d]], 30, $Aborted]];
  Do[
    Module[{dv = TimeConstrained[Quiet[N[d /. x -> p, 40]], 10, $Aborted],
            fv = Quiet[N[f /. x -> p, 40]], ok},
      ok = If[NumberQ[dv] && NumberQ[fv],
        If[Abs[dv] <= 10^-TOL Max[1, Abs[fv]], "ok", "bad"], "na"];
      Print["@@NUM ", InputForm[p], " ", ok, " ",
        InputForm[If[NumberQ[dv], N[Abs[dv], 6], dv]]]],
    {p, pts}]];
Print["@@DONE"];
"""


def run_mathilda(mbin, case, method, timeout, tmpdir):
    body = (RUN_TMPL.replace("INTEGRAND", case["f"])
            .replace("POINTS", "{" + ", ".join(case["points"]) + "}")
            .replace("WORKER", "Integrate`" + method)
            .replace("TIMEOUT", str(timeout)).replace("TOL", str(TOL_EXP)))
    path = os.path.join(tmpdir, f"{case['id']}_{method}.m")
    mathilda_script(path, body)
    out, err, rc, wall, to = run_proc([mbin, "-file", path],
                                      timeout + VERIFY_BUDGET)
    rec = {"wall": round(wall, 3)}

    def grab(tag):
        m = re.search(rf"^@@{tag} (.*)$", out, re.M)
        return m.group(1).strip() if m else None
    t = grab("TIME")
    rec["time"] = float(re.sub(r"`.*", "", t)) if t and re.match(r"[\d.]", t) else None
    rec["result"] = grab("RES")
    st = grab("STATUS")
    if to and st is None:
        rec.update(status="TIMEOUT", detail="process timeout")
        rec["time"] = rec["time"] or round(wall, 3)
        return rec
    if st is None:
        rec.update(status="ERROR",
                   detail=f"rc={rc}; " + (err.strip()[-400:] or out.strip()[-400:]))
        return rec
    if st == "timeout":
        rec.update(status="TIMEOUT", detail="TimeConstrained")
    elif st == "not elementary":
        rec.update(status="NONELEM", detail=rec["result"])
    elif st in ("failed", "unevaluated"):
        rec.update(status="DECLINE", detail=rec["result"])
    elif st == "elem":
        pzq = grab("PZQ")
        nums = re.findall(r"^@@NUM (\S+) (ok|bad|na) (.*)$", out, re.M)
        rec["pzq"] = pzq
        rec["numeric"] = [list(n) for n in nums]
        # The 40-digit numeric check is decisive: a PossibleZeroQ False
        # contradicted by three exact-point checks is a PossibleZeroQ false
        # negative (reported separately), not a wrong antiderivative.
        bad = any(n[1] == "bad" for n in nums)
        okn = sum(n[1] == "ok" for n in nums)
        if bad:
            rec["status"] = "WRONG"
        elif okn == len(nums) and okn >= 2:
            rec["status"] = "ELEM"
            if pzq == "False":
                rec["pzq_false_negative"] = True
                rec["detail"] = "PossibleZeroQ[D[F, x] - f] is False, numeric check passes"
        elif pzq == "True":
            rec["status"] = "ELEM"
            rec["detail"] = "verified by PossibleZeroQ only (numeric check incomplete)"
        else:
            rec["status"] = "WRONG"
            rec["detail"] = ("could not be verified (PossibleZeroQ %s, %d/%d points)"
                             % (pzq, okn, len(nums)))
    else:
        rec.update(status="ERROR", detail=f"unknown status {st}")
    return rec


# --------------------------------------------------------------------------
# FriCAS
# --------------------------------------------------------------------------

def fricas_script(case, outfile):
    f = to_fricas(case["f"])
    L = [")set output algebra off", ")set messages type off",
         f'tf : TextFile := open("{outfile}", "output")',
         'writeLine!(tf, "START")', "flush(tf)",
         f"f := {f}",
         ")set messages time on",
         "r := integrate(f, x)",
         ")set messages time off",
         'writeLine!(tf, concat("RES ", unparse(r::InputForm)))', "flush(tf)"]
    for p in case["points"]:
        L.append(f"c := complexNumeric(eval(D(r, x) - f, x, {p}), 40)")
        L.append(f"cf := complexNumeric(eval(f, x, {p}), 40)")
        L.append(f'writeLine!(tf, concat(["NUM {p} ", '
                 f'unparse(convert(sqrt(norm(c)))@InputForm), " ", '
                 f'unparse(convert(sqrt(norm(cf)))@InputForm)]))')
        L.append("flush(tf)")
    L += ['writeLine!(tf, "DONE")', "close!(tf)", ")quit"]
    return "\n".join(L) + "\n"


def fricas_float(s):
    s = s.strip()
    m = re.fullmatch(r"float\((-?\d+),\s*(-?\d+),\s*(\d+)\)", s)
    if m:
        mant, ex, base = int(m.group(1)), int(m.group(2)), int(m.group(3))
        return Fraction(mant) * Fraction(base) ** ex
    try:
        return Fraction(s)
    except (ValueError, ZeroDivisionError):
        return None


def run_fricas(fbin, case, timeout, tmpdir):
    outfile = os.path.join(tmpdir, f"{case['id']}_fricas.txt")
    if os.path.exists(outfile):
        os.unlink(outfile)
    out, err, rc, wall, to = run_proc([fbin, "-nosman"], timeout,
                                      stdin=fricas_script(case, outfile))
    rec = {"wall": round(wall, 3)}
    m = re.search(r"Time: .*?= ([\d.]+) sec", out)
    rec["time"] = float(m.group(1)) if m else None
    lines = open(outfile).read().splitlines() if os.path.exists(outfile) else []
    res = [ln[4:] for ln in lines if ln.startswith("RES ")]
    rec["result"] = res[0] if res else None
    if not res:
        if to:
            rec.update(status="TIMEOUT", time=rec["time"] or round(wall, 3))
        else:
            em = re.search(r">> (?:Error detected within library code:\s*)?(.*?)\n",
                           out)
            rec.update(status="ERROR",
                       detail=(em.group(1).strip() if em else out.strip()[-300:]))
        return rec
    r = res[0]
    if r == "r":            # integrate signalled an error, r was never bound
        em = re.search(r">> (?:Error detected within library code:\s*)?(.*?)\n\s*\n",
                       out, re.S)
        rec.update(status="ERROR", result=None,
                   detail=(re.sub(r"\s+", " ", em.group(1)).strip() if em
                           else "integrate failed"))
        return rec
    if "integral(" in r:
        rec["status"] = "NONELEM"
        return rec
    if r.startswith("["):
        rec.update(status="ERROR", detail="list of results (parameter case split)")
        return rec
    nums = []
    for ln in lines:
        mm = re.match(r"NUM (\S+) (float\([^)]*\)|\S+) (float\([^)]*\)|\S+)$", ln)
        if mm:
            dv, fv = fricas_float(mm.group(2)), fricas_float(mm.group(3))
            if dv is None or fv is None:
                nums.append([mm.group(1), "na", mm.group(2)])
            else:
                ok = dv <= Fraction(1, 10 ** TOL_EXP) * max(1, fv)
                nums.append([mm.group(1), "ok" if ok else "bad",
                             f"{float(dv):.3g}"])
    rec["numeric"] = nums
    if any(n[1] == "bad" for n in nums):
        rec["status"] = "BAD"
    elif any(n[1] == "ok" for n in nums):
        rec["status"] = "ELEM"
    else:
        rec.update(status="ELEM", detail="unverified (numeric check failed to run)")
    return rec


# --------------------------------------------------------------------------
# Classification and reports
# --------------------------------------------------------------------------

def classify(m, fr, expect):
    ms, fs = m["status"], fr["status"]
    f_elem = fs == "ELEM"
    if ms == "WRONG":
        return "WRONG"
    if ms == "NONELEM":
        if f_elem or expect == "elem":
            return "DISAGREE"
        return "AGREE-NONELEM" if fs == "NONELEM" else "NONELEM-UNCONFIRMED"
    if ms == "ELEM":
        return "AGREE-ELEM" if f_elem else "MATHILDA-ONLY"
    if ms in ("DECLINE", "TIMEOUT") and f_elem:
        return "FRICAS-ONLY"
    if ms == "DECLINE":
        return "NEITHER"
    return ms            # TIMEOUT / ERROR


def trunc(s, n=160):
    s = "" if s is None else str(s)
    return s if len(s) <= n else s[:n] + "..."


def md_cell(s, n=120):
    return trunc(s, n).replace("|", "\\|").replace("\n", " ")


def tstats(xs):
    xs = [v for v in xs if v is not None]
    if not xs:
        return "-", "-", "-", 0
    return (f"{statistics.median(xs):.3f}", f"{statistics.mean(xs):.3f}",
            f"{max(xs):.3f}", len(xs))


def write_reports(recs, args, elapsed):
    counts = {mt: {c: 0 for c in CLASSES} for mt in METHODS}
    fams = {}
    for r in recs:
        for mt in METHODS:
            counts[mt][r["class"][mt]] += 1
            fams.setdefault(r["family"], {m2: {} for m2 in METHODS})
            d = fams[r["family"]][mt]
            d[r["class"][mt]] = d.get(r["class"][mt], 0) + 1
    L = ["# Trager vs FriCAS stress report", "",
         f"Generated {time.strftime('%Y-%m-%d %H:%M')} by "
         f"`tools/trager_fricas_stress.py` over `tests/data/trager_stress_corpus.m` "
         f"({len(recs)} cases, timeout {args.timeout} s per case and system, "
         f"{args.jobs} jobs, {elapsed:.0f} s wall).", "",
         "Classes: see the docstring of the tool.  WRONG and DISAGREE are hard "
         "failures; FRICAS-ONLY ratchets against `BASELINE`.", "",
         "## Summary", "",
         "| class | " + " | ".join(METHODS) + " |",
         "|---|" + "---:|" * len(METHODS)]
    for c in CLASSES:
        L.append(f"| {c} | " + " | ".join(str(counts[mt][c]) for mt in METHODS) + " |")
    fst = {}
    for r in recs:
        fst[r["fricas"]["status"]] = fst.get(r["fricas"]["status"], 0) + 1
    L += ["", "FriCAS outcomes: " + ", ".join(f"{k} {v}" for k, v in sorted(fst.items()))
          + " (BAD = FriCAS antiderivative failed its own derivative check).", ""]
    L += ["## Per family", "",
          "| family | n | Trager | TragerKauers |", "|---|---:|---|---|"]
    for fam in sorted(fams):
        n = sum(fams[fam]["Trager"].values())
        cells = [", ".join(f"{k} {v}" for k, v in sorted(fams[fam][mt].items()))
                 for mt in METHODS]
        L.append(f"| {fam} | {n} | {cells[0]} | {cells[1]} |")
    # timing
    L += ["", "## Timing (seconds; Mathilda: AbsoluteTiming of the worker; "
          "FriCAS: its own `Time:` total for integrate)", "",
          "| system | median | mean | max | n |", "|---|---:|---:|---:|---:|"]
    for mt in METHODS:
        L.append(f"| Mathilda {mt} | " + " | ".join(
            map(str, tstats([r["mathilda"][mt].get("time") for r in recs]))) + " |")
    L.append("| FriCAS | " + " | ".join(
        map(str, tstats([r["fricas"].get("time") for r in recs]))) + " |")
    slow = []
    for r in recs:
        ft = r["fricas"].get("time")
        for mt in METHODS:
            t = r["mathilda"][mt].get("time")
            if t is not None and ft is not None and t > 0.1 and \
                    t > 10 * max(ft, 0.01):
                slow.append((r["id"], mt, t, ft))
    L += ["", f"### Mathilda > 10x slower than FriCAS ({len(slow)})", ""]
    if slow:
        L += ["| id | method | Mathilda s | FriCAS s |", "|---|---|---:|---:|"]
        L += [f"| {a} | {b} | {c:.3f} | {d:.2f} |" for a, b, c, d in slow]
    # problem lists
    L += ["", "## Hard failures and gaps", ""]
    for title, pred in (("WRONG", lambda r, mt: r["class"][mt] == "WRONG"),
                        ("DISAGREE", lambda r, mt: r["class"][mt] == "DISAGREE"),
                        ("FRICAS-ONLY", lambda r, mt: r["class"][mt] == "FRICAS-ONLY")):
        rows = [(r, mt) for r in recs for mt in METHODS if pred(r, mt)]
        L += [f"### {title} ({len(rows)})", ""]
        for r, mt in rows:
            m = r["mathilda"][mt]
            L.append(f"- `{r['id']}` [{mt}] `{r['f']}` -- Mathilda: {m['status']} "
                     f"`{md_cell(m.get('result'), 200)}`; FriCAS: "
                     f"{r['fricas']['status']} `{md_cell(r['fricas'].get('result'), 200)}`")
        L.append("")
    pzf = [(r, mt) for r in recs for mt in METHODS
           if r["mathilda"][mt].get("pzq_false_negative")]
    L += [f"### PossibleZeroQ false negatives ({len(pzf)})", "",
          "Verified antiderivatives (40-digit check at every point) for which "
          "`PossibleZeroQ[D[F, x] - f]` returned False -- a PossibleZeroQ bug, "
          "not an integrator bug.", ""]
    for r, mt in pzf:
        L.append(f"- `{r['id']}` [{mt}] `{r['f']}`: "
                 f"`{md_cell(r['mathilda'][mt].get('result'), 200)}`")
    L.append("")
    badf = [r for r in recs if r["fricas"]["status"] == "BAD"]
    L += [f"### FriCAS answers failing their derivative check ({len(badf)})", ""]
    for r in badf:
        L.append(f"- `{r['id']}` `{r['f']}`: `{md_cell(r['fricas'].get('result'), 200)}` "
                 f"numeric {r['fricas'].get('numeric')}")
    L += ["", "## All cases", "",
          "| id | family | integrand | expect | Trager | s | TragerKauers | s | "
          "FriCAS | s | class T | class TK |",
          "|---|---|---|---|---|---:|---|---:|---|---:|---|---|"]
    for r in recs:
        mT, mK, fr = r["mathilda"]["Trager"], r["mathilda"]["TragerKauers"], r["fricas"]

        def ts(v):
            return "-" if v is None else f"{v:.3f}"
        L.append(f"| {r['id']} | {r['family']} | `{md_cell(r['f'], 90)}` | "
                 f"{r['expect']} | {mT['status']} | {ts(mT.get('time'))} | "
                 f"{mK['status']} | {ts(mK.get('time'))} | {fr['status']} | "
                 f"{ts(fr.get('time'))} | {r['class']['Trager']} | "
                 f"{r['class']['TragerKauers']} |")
    with open(args.report_md, "w") as fh:
        fh.write("\n".join(L) + "\n")
    with open(args.report_json, "w") as fh:
        json.dump({"generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
                   "timeout": args.timeout, "counts": counts,
                   "baseline": BASELINE, "cases": recs}, fh, indent=1)
    return counts


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--filter", default=None, help="regex on id, family or integrand")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--timeout", type=int, default=60)
    ap.add_argument("--corpus", default=CORPUS)
    ap.add_argument("--report-md", default=None)
    ap.add_argument("--report-json", default=None)
    ap.add_argument("--emit-random", action="store_true",
                    help="print the seeded corpus sections and exit")
    args = ap.parse_args()

    if args.emit_random:
        print(",\n".join(corpus_line(*t) for t in gen_seeded()))
        return 0

    fbin = find_fricas()
    if not fbin:
        print("trager_fricas_stress: fricas not found -- skipping (exit 0)")
        return 0
    mbin = find_mathilda()
    if not os.path.exists(mbin):
        print(f"trager_fricas_stress: Mathilda binary {mbin} not found")
        return 1

    cases = load_corpus(args.corpus)
    if args.filter:
        rx = re.compile(args.filter)
        cases = [c for c in cases if rx.search(c["id"]) or rx.search(c["family"])
                 or rx.search(c["f"])]
    if args.limit:
        cases = cases[:args.limit]
    partial = bool(args.filter or args.limit)
    # A partial run must not overwrite the checked-in full report.
    if args.report_md is None:
        args.report_md = REPORT_MD if not partial else \
            os.path.join(tempfile.gettempdir(), "trager_stress_report.md")
    if args.report_json is None:
        args.report_json = REPORT_JSON if not partial else \
            os.path.join(tempfile.gettempdir(), "trager_stress_report.json")

    t0 = time.time()
    tmpdir = tempfile.mkdtemp(prefix="trager_stress_")
    print(f"trager_fricas_stress: {len(cases)} cases, {args.jobs} jobs, "
          f"timeout {args.timeout} s; scratch {tmpdir}", flush=True)
    mathilda_points(mbin, cases, tmpdir)

    lock = threading.Lock()
    done = [0]

    def work(c):
        rec = {"id": c["id"], "family": c["family"], "f": c["f"],
               "expect": c["expect"], "points": c["points"],
               "fricas_input": to_fricas(c["f"]), "mathilda": {}}
        for mt in METHODS:
            rec["mathilda"][mt] = run_mathilda(mbin, c, mt, args.timeout, tmpdir)
        rec["fricas"] = run_fricas(fbin, c, args.timeout, tmpdir)
        rec["class"] = {mt: classify(rec["mathilda"][mt], rec["fricas"], c["expect"])
                        for mt in METHODS}
        with lock:
            done[0] += 1
            flag = " <<<" if any(rec["class"][mt] in HARD + ("FRICAS-ONLY",)
                                 for mt in METHODS) else ""
            print(f"[{done[0]:3d}/{len(cases)}] {c['id']:<22} "
                  f"T {rec['class']['Trager']:<20} TK {rec['class']['TragerKauers']:<20}"
                  f" F {rec['fricas']['status']}{flag}", flush=True)
        return rec

    with cf.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        recs = list(ex.map(work, cases))
    shutil.rmtree(tmpdir, ignore_errors=True)
    counts = write_reports(recs, args, time.time() - t0)

    print("\n%-20s %8s %13s" % ("class", *METHODS))
    for c in CLASSES:
        print("%-20s %8d %13d" % (c, counts["Trager"][c], counts["TragerKauers"][c]))
    print(f"reports: {args.report_md}, {args.report_json}")

    fail = False
    for mt in METHODS:
        for c in HARD:
            if counts[mt][c]:
                fail = True
                ids = [r["id"] for r in recs if r["class"][mt] == c]
                print(f"FAIL [{mt}] {c}: {', '.join(ids)}")
        n = counts[mt]["FRICAS-ONLY"]
        if not partial and n > BASELINE[mt]:
            fail = True
            print(f"FAIL [{mt}] FRICAS-ONLY {n} > BASELINE {BASELINE[mt]}")
        elif not partial and n < BASELINE[mt]:
            print(f"note [{mt}] FRICAS-ONLY {n} < BASELINE {BASELINE[mt]}: "
                  f"lower the BASELINE (ratchet)")
    if args.check:
        return 1 if fail else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
