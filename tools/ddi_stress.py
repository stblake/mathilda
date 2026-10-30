#!/usr/bin/env python3
"""
ddi_stress.py -- stress test of Mathilda's algebraic integrators

    Integrate[f, x, Method -> "Trager"]              (Integrate`Trager)
    Integrate[f, x, Method -> "ParallelMixedTower"]  (Integrate`ParallelMixedTower)

on the integrals of S. Blake, "Diabolically Difficult Integrals", Volume 1,
transcribed in tests/data/ddi_vol1_corpus.m, with FriCAS as a third opinion.

It is a sibling of tools/trager_fricas_stress.py and reuses its process
runner, FriCAS float parser and report helpers.  Differences:

  * the corpus syntax is wider: symbolic parameters (single lowercase
    letters), nested and multiple radicals, Log/Exp/Sin/Cos/Tan/ArcTan;
  * every record carries a book reference;
  * antiderivatives of parametric integrands are verified numerically after
    substituting fixed generic rational values for the parameters (PARAMS),
    at three points x where the integrand is real and finite;
  * FriCAS answers that split on the sign of a parameter (a list) are
    verified branch by branch -- one verified branch counts as an answer.

Per method each case is classified as in trager_fricas_stress.py
(AGREE-ELEM, AGREE-NONELEM, MATHILDA-ONLY, FRICAS-ONLY, DISAGREE, WRONG,
NONELEM-UNCONFIRMED, NEITHER, TIMEOUT, ERROR); in addition a verified Mathilda
antiderivative for an integrand the book calls non-elementary is flagged
BOOK-CONFLICT (the verified antiderivative is itself a proof).  An answer whose
derivative is +f on part of the real domain and -f on the rest (a radical read
on the wrong branch there) is BRANCH; one whose check could not finish within
the budget is UNVERIFIED.  A FriCAS answer in special functions (erf, dilog,
...) counts as FriCAS deciding the integral is not elementary.

Usage (from the repository root):
    python3 tools/ddi_stress.py                       # full run
    python3 tools/ddi_stress.py --filter nested --jobs 8 --timeout 60
    python3 tools/ddi_stress.py --no-fricas
    python3 tools/ddi_stress.py --merge --no-fricas --methods ParallelMixedTower
                                    # re-run one method, keep the other columns
    python3 tools/ddi_stress.py --reclassify   # after editing expect values
Writes tests/data/ddi_vol1_report.json.  Python 3 standard library only.
"""
import argparse
import concurrent.futures as cf
import json
import os
import re
import sys
import tempfile
import threading
import time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trager_fricas_stress as tfs  # noqa: E402  (runner + helpers)

REPO = tfs.REPO
CORPUS = os.path.join(REPO, "tests", "data", "ddi_vol1_corpus.m")
REPORT_JSON = os.path.join(REPO, "tests", "data", "ddi_vol1_report.json")

METHODS = ("Trager", "ParallelMixedTower")
TOL_EXP = 20
VERIFY_BUDGET = 90

# Generic values substituted for parameters in the numeric checks only.
PARAMS = {"a": "3/2", "b": "5/3", "c": "7/4", "d": "9/5", "e": "11/6",
          "k": "4/3", "m": "5/2", "n": "7/3", "p": "8/5", "q": "10/7",
          "r": "12/7", "s": "13/8", "t": "14/9"}

FUNCS = ("Sqrt", "Log", "Exp", "Sin", "Cos", "Tan", "ArcTan")
FRICAS_FUNCS = {"Sqrt": "sqrt", "Log": "log", "Exp": "exp", "Sin": "sin",
                "Cos": "cos", "Tan": "tan", "ArcTan": "atan"}

LINE_RE = re.compile(r'^\s*\{\s*"([^"]+)"\s*,\s*"([^"]+)"\s*,\s*"([^"]+)"\s*,'
                     r'\s*"([^"]*)"\s*,\s*"([^"]*)"\s*\}\s*,?\s*$')
TOKEN_RE = re.compile(r"\s*(" + "|".join(f + r"\[" for f in FUNCS) +
                      r"|\d+|[a-z]|[-+*/^()\]])")


def tokens(f):
    """Tokenise the corpus syntax; raise ValueError on anything else (in
    particular an implicit product, which FriCAS would read as a call)."""
    pos, prev, depth, out = 0, None, [], []
    for m in TOKEN_RE.finditer(f):
        if m.start() != pos:
            raise ValueError(f"bad character at {pos}: {f!r}")
        tok = m.group(1)
        operand = tok.isdigit() or tok.isalpha() or tok == "(" or tok.endswith("[")
        if operand and prev is not None and (prev in (")", "]") or prev.isdigit()
                                             or (prev.isalpha() and len(prev) == 1)):
            raise ValueError(f"implicit product before {tok!r}: {f!r}")
        if tok == "(" or tok.endswith("["):
            depth.append(tok)
        elif tok in (")", "]"):
            if not depth or (depth.pop() == "(") != (tok == ")"):
                raise ValueError(f"unbalanced {tok!r}: {f!r}")
        if len(tok) == 1 and tok.isalpha() and tok != "x" and tok not in PARAMS:
            raise ValueError(f"unknown parameter {tok!r}: {f!r}")
        out.append(tok)
        prev = tok
        pos = m.end()
    if pos != len(f.rstrip()) or depth:
        raise ValueError(f"unparsable or unbalanced: {f!r}")
    if "x" not in out:
        raise ValueError(f"no x: {f!r}")
    return out


def params_of(f):
    return sorted({t for t in tokens(f) if len(t) == 1 and t.isalpha() and t != "x"})


def to_fricas(f):
    out = []
    for t in tokens(f):
        out.append(FRICAS_FUNCS[t[:-1]] + "(" if t.endswith("[") else
                   ")" if t == "]" else t)
    return "".join(out)


def load_corpus(path):
    cases, seen = [], set()
    text = open(path).read()
    text = re.sub(r"\(\*.*?\*\)", lambda m: "\n" * m.group(0).count("\n"),
                  text, flags=re.S)
    for ln, line in enumerate(text.split("\n"), 1):
        m = LINE_RE.match(line)
        if not m:
            continue
        cid, fam, f, expect, ref = m.groups()
        if cid in seen:
            raise SystemExit(f"{path}:{ln}: duplicate id {cid}")
        seen.add(cid)
        try:
            ps = params_of(f)
        except ValueError as e:
            raise SystemExit(f"{path}:{ln}: {e}")
        cases.append({"id": cid, "family": fam, "f": f, "expect": expect,
                      "ref": ref, "params": ps})
    return cases


def subst_rules(case):
    """Mathilda replacement list for the parameters, e.g. {a -> 3/2}."""
    return "{" + ", ".join(f"{p} -> {PARAMS[p]}" for p in case["params"]) + "}"


# --------------------------------------------------------------------------
# Mathilda
# --------------------------------------------------------------------------

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
    """Three evaluation points per case (parameters at PARAMS), batched."""
    pre = POINTS_TMPL.replace("CANDS", "{" + ", ".join(tfs.CANDIDATES) + "}")

    def batch(cs):
        body = pre + "".join(
            f'Print["@@PTS {c["id"]} ", InputForm[TimeConstrained['
            f'Stress`Pts[({c["f"]}) /. {subst_rules(c)}], 20, {{}}]]];\n'
            for c in cs)
        path = os.path.join(tmpdir, f"pts_{cs[0]['id']}.m")
        tfs.mathilda_script(path, body)
        out, _, _, _, _ = tfs.run_proc([mbin, "-file", path], 30 + 20 * len(cs))
        res = {}
        for m in re.finditer(r"^@@PTS (\S+) \{(.*)\}\s*$", out, re.M):
            res[m.group(1)] = [s.strip() for s in m.group(2).split(",") if s.strip()]
        return res

    res = {}
    for i in range(0, len(cases), 40):
        res.update(batch(cases[i:i + 40]))
    for c in cases:
        if c["id"] not in res:
            res.update(batch([c]))
        c["points"] = res.get(c["id"]) or tfs.CANDIDATES[:3]


RUN_TMPL = r"""
f = (INTEGRAND);
pts = POINTS;
prm = PRM;
{t, r} = AbsoluteTiming[TimeConstrained[Quiet[WORKER[f, x]], TIMEOUT, $Aborted]];
Print["@@TIME ", InputForm[N[t]]];
Print["@@RES ", InputForm[r]];
Print["@@LEAVES ", LeafCount[r]];
status = Which[
  r === $Aborted, "timeout",
  ListQ[r] && Length[r] >= 2 && StringQ[r[[1]]], r[[1]],
  !FreeQ[r, Integrate | WORKER], "unevaluated",
  True, "elem"];
Print["@@STATUS ", status];
If[status === "elem",
  (* numeric check first, on the antiderivative with the parameters already
     instantiated (differentiating the symbolic answer and substituting later
     can exhaust the process budget on large parametric answers); "neg" marks
     D[F] = -f exactly, the signature of a radical read on the wrong branch *)
  dn = TimeConstrained[Quiet[D[r /. prm, x] - (f /. prm)], 30, $Aborted];
  Do[
    Module[{dv = TimeConstrained[Quiet[N[dn /. x -> p, 40]], 15, $Aborted],
            fv = Quiet[N[(f /. prm) /. x -> p, 40]], ok},
      ok = If[NumberQ[dv] && NumberQ[fv],
        Which[Abs[dv] <= 10^-TOL Max[1, Abs[fv]], "ok",
              Abs[dv + 2 fv] <= 10^-TOL Max[1, Abs[fv]], "neg", True, "bad"], "na"];
      Print["@@NUM ", InputForm[p], " ", ok, " ",
        InputForm[If[NumberQ[dv], N[Abs[dv], 6], dv]]]],
    {p, pts}];
  Print["@@PZQ ", TimeConstrained[Quiet[PossibleZeroQ[D[r, x] - f]], 30, $Aborted]]];
Print["@@DONE"];
"""


def run_mathilda(mbin, case, method, timeout, tmpdir):
    body = (RUN_TMPL.replace("INTEGRAND", case["f"])
            .replace("POINTS", "{" + ", ".join(case["points"]) + "}")
            .replace("PRM", subst_rules(case))
            .replace("WORKER", "Integrate`" + method)
            .replace("TIMEOUT", str(timeout)).replace("TOL", str(TOL_EXP)))
    path = os.path.join(tmpdir, f"{case['id']}_{method}.m")
    tfs.mathilda_script(path, body)
    out, err, rc, wall, to = tfs.run_proc([mbin, "-file", path],
                                          timeout + VERIFY_BUDGET)
    rec = {"wall": round(wall, 3)}

    def grab(tag):
        m = re.search(rf"^@@{tag} (.*)$", out, re.M)
        return m.group(1).strip() if m else None
    t = grab("TIME")
    rec["time"] = float(re.sub(r"`.*", "", t)) if t and re.match(r"[\d.]", t) else None
    rec["result"] = grab("RES")
    lv = grab("LEAVES")
    rec["leaves"] = int(lv) if lv and lv.isdigit() else None
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
    elif st == "elem":
        pzq = grab("PZQ")
        nums = re.findall(r"^@@NUM (\S+) (ok|bad|neg|na) (.*)$", out, re.M)
        rec["pzq"] = pzq
        rec["numeric"] = [list(n) for n in nums]
        bad = any(n[1] == "bad" for n in nums)
        okn = sum(n[1] == "ok" for n in nums)
        negn = sum(n[1] == "neg" for n in nums)
        if bad:
            rec["status"] = "WRONG"
        elif negn and okn:
            # D[F] = +f on part of the real domain and -f on the rest: a
            # principal-branch answer that is wrong on some real interval
            rec["status"] = "BRANCH"
            rec["detail"] = "D[F] = f at %d point(s), -f at %d" % (okn, negn)
        elif negn:
            rec["status"] = "WRONG"
            rec["detail"] = "D[F] = -f at every sample (global sign)"
        elif okn == len(nums) and okn >= 2:
            rec["status"] = "ELEM"
            if pzq == "False":
                rec["detail"] = "PossibleZeroQ[D[F, x] - f] is False, numeric check passes"
        elif pzq == "True":
            rec["status"] = "ELEM"
            rec["detail"] = "verified by PossibleZeroQ only (numeric check incomplete)"
        else:
            rec["status"] = "UNVERIFIED"
            rec["detail"] = ("could not be verified (PossibleZeroQ %s, %d/%d points)"
                             % (pzq, okn, len(nums)))
    else:   # "failed", "unevaluated", or any other worker verdict string
        rec.update(status="DECLINE", detail=rec["result"])
    return rec


# --------------------------------------------------------------------------
# FriCAS
# --------------------------------------------------------------------------

def fricas_script(case, outfile):
    f = to_fricas(case["f"])
    sub = ", ".join(f"{p} = {PARAMS[p]}" for p in case["params"])

    def at(p):
        return f"[x = {p}" + (", " + sub if sub else "") + "]"
    L = [")set output algebra off", ")set messages type off",
         f'tf : TextFile := open("{outfile}", "output")',
         'writeLine!(tf, "START")', "flush(tf)",
         f"f := {f}",
         ")set messages time on",
         "r := integrate(f, x)",
         ")set messages time off",
         'writeLine!(tf, concat("RES ", unparse(r::InputForm)))', "flush(tf)",
         # a parameter case split is a List(Expression Integer): check each
         "rl := (if r case List(Expression(Integer)) then r::List(Expression(Integer))"
         " else [r::Expression(Integer)])"]
    for p in case["points"]:
        L.append(f"for g in rl repeat (c := complexNumeric(eval(D(g, x) - f, {at(p)}), 40); "
                 f"cf := complexNumeric(eval(f, {at(p)}), 40); "
                 f'writeLine!(tf, concat(["NUM {p} ", '
                 f'unparse(convert(sqrt(norm(c)))@InputForm), " ", '
                 f'unparse(convert(sqrt(norm(cf)))@InputForm)])); flush(tf))')
    L += ['writeLine!(tf, "DONE")', "close!(tf)", ")quit"]
    return "\n".join(L) + "\n"


FRICAS_SPECIAL_RE = re.compile(r"\b(erfi?|Ei|li|Si|Ci|dilog|polylog|fresnel[SC]|"
                               r"elliptic[A-Za-z]*|lambertW|Gamma|airy[A-Za-z]*|"
                               r"bessel[A-Za-z]*)\(")


def run_fricas(fbin, case, timeout, tmpdir):
    outfile = os.path.join(tmpdir, f"{case['id']}_fricas.txt")
    if os.path.exists(outfile):
        os.unlink(outfile)
    out, err, rc, wall, to = tfs.run_proc([fbin, "-nosman"], timeout,
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
            em = re.search(r">> (?:Error detected within library code:\s*)?(.*?)\n", out)
            rec.update(status="ERROR",
                       detail=(em.group(1).strip() if em else out.strip()[-300:]))
        return rec
    r = res[0]
    if r == "r":
        em = re.search(r">> (?:Error detected within library code:\s*)?(.*?)\n\s*\n",
                       out, re.S)
        rec.update(status="ERROR", result=None,
                   detail=(re.sub(r"\s+", " ", em.group(1)).strip() if em
                           else "integrate failed"))
        return rec
    if "integral(" in r:
        rec["status"] = "NONELEM"
        return rec
    sf = FRICAS_SPECIAL_RE.findall(r)
    if sf:
        # FriCAS left the elementary functions (erf, dilog, Ei, ...): its
        # integrate decided the integral is NOT elementary, and a verified
        # special-function answer must not count as an elementary one
        rec.update(status="NONELEM", detail="special function: " + ", ".join(sorted(set(sf))))
        return rec
    # NUM lines: one per (point, branch).  A branch is good when all its
    # points pass; with a single branch this is the usual check.
    per_point = {}
    for ln in lines:
        mm = re.match(r"NUM (\S+) (float\([^)]*\)|\S+) (float\([^)]*\)|\S+)$", ln)
        if mm:
            dv, fv = tfs.fricas_float(mm.group(2)), tfs.fricas_float(mm.group(3))
            st = ("na" if dv is None or fv is None else
                  "ok" if dv <= Fraction(1, 10 ** TOL_EXP) * max(1, fv) else "bad")
            per_point.setdefault(mm.group(1), []).append(st)
    rec["numeric"] = per_point
    nbr = max((len(v) for v in per_point.values()), default=0)
    branch_ok = [all(v[i] == "ok" for v in per_point.values() if len(v) > i)
                 and any(len(v) > i for v in per_point.values())
                 for i in range(nbr)]
    if r.startswith("["):
        rec["split"] = True
    if any(branch_ok):
        rec["status"] = "ELEM"
    elif any("bad" in v for v in per_point.values()):
        rec["status"] = "BAD"
    else:
        rec.update(status="ELEM", detail="unverified (numeric check failed to run)")
    return rec


# --------------------------------------------------------------------------
# Classification
# --------------------------------------------------------------------------

CLASSES = ("AGREE-ELEM", "AGREE-NONELEM", "MATHILDA-ONLY", "FRICAS-ONLY",
           "DISAGREE", "WRONG", "BOOK-CONFLICT", "NONELEM-UNCONFIRMED",
           "NEITHER", "TIMEOUT", "ERROR")


def classify(m, fr, expect):
    ms, fs = m["status"], fr["status"]
    f_elem = fs == "ELEM"
    if ms in ("WRONG", "BRANCH", "UNVERIFIED"):
        return ms
    if ms == "NONELEM":
        if f_elem or expect == "elem":
            return "DISAGREE"
        return "AGREE-NONELEM" if fs == "NONELEM" or expect == "nonelem" \
            else "NONELEM-UNCONFIRMED"
    if ms == "ELEM":
        if expect == "nonelem":
            return "BOOK-CONFLICT"
        return "AGREE-ELEM" if f_elem else "MATHILDA-ONLY"
    if ms in ("DECLINE", "TIMEOUT") and f_elem:
        return "FRICAS-ONLY"
    if ms == "DECLINE":
        return "NEITHER"
    return ms


def reclassify(corpus, out):
    d = json.load(open(out))
    expect = {c["id"]: c["expect"] for c in load_corpus(corpus)}
    for c in d["cases"]:
        c["expect"] = expect.get(c["id"], c["expect"])
        fr = c["fricas"]
        sf = FRICAS_SPECIAL_RE.findall(fr.get("result") or "")
        if fr["status"] in ("ELEM", "BAD") and sf:
            fr.update(status="NONELEM",
                      detail="special function: " + ", ".join(sorted(set(sf))))
        for m in c["methods"].values():
            m["class"] = classify(m, fr, c["expect"])
    with open(out, "w") as fh:
        json.dump(d, fh, indent=1)
    print(f"reclassified {len(d['cases'])} cases in {out}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus", default=CORPUS)
    ap.add_argument("--out", default=REPORT_JSON)
    ap.add_argument("--filter", default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--jobs", type=int, default=max(2, (os.cpu_count() or 4) - 2))
    ap.add_argument("--timeout", type=int, default=60)
    ap.add_argument("--fricas-timeout", type=int, default=60)
    ap.add_argument("--no-fricas", action="store_true")
    ap.add_argument("--methods", default=",".join(METHODS))
    ap.add_argument("--merge", action="store_true",
                    help="re-run only --methods (and FriCAS unless --no-fricas) on the "
                         "selected cases and splice the results into the existing --out "
                         "report, keeping every other column and case")
    ap.add_argument("--reclassify", action="store_true",
                    help="run nothing: re-read expect values from --corpus, re-apply the "
                         "FriCAS special-function rule to the stored answers, and "
                         "reclassify every case of --out in place")
    args = ap.parse_args()

    if args.reclassify:
        reclassify(args.corpus, args.out)
        return
    cases = load_corpus(args.corpus)
    if args.filter:
        cases = [c for c in cases if re.search(args.filter, c["id"] + " " + c["family"])]
    if args.limit:
        cases = cases[:args.limit]
    methods = [m for m in args.methods.split(",") if m]
    mbin = tfs.find_mathilda()
    fbin = None if args.no_fricas else tfs.find_fricas()
    print(f"{len(cases)} cases, methods {methods}, fricas {fbin}, jobs {args.jobs}",
          file=sys.stderr)

    t0 = time.time()
    tmpdir = tempfile.mkdtemp(prefix="ddi_stress_")
    mathilda_points(mbin, cases, tmpdir)

    jobs = [(c, mt) for c in cases for mt in methods]
    if fbin:
        jobs += [(c, "fricas") for c in cases]
    lock, done = threading.Lock(), [0]

    def work(job):
        c, mt = job
        if mt == "fricas":
            r = run_fricas(fbin, c, args.fricas_timeout, tmpdir)
        else:
            r = run_mathilda(mbin, c, mt, args.timeout, tmpdir)
        with lock:
            done[0] += 1
            print(f"[{done[0]}/{len(jobs)}] {c['id']} {mt}: {r['status']} "
                  f"{r.get('time')}", file=sys.stderr, flush=True)
        return c["id"], mt, r

    results = {}
    with cf.ThreadPoolExecutor(args.jobs) as ex:
        for cid, mt, r in ex.map(work, jobs):
            results.setdefault(cid, {})[mt] = r

    old = json.load(open(args.out)) if args.merge else None
    prev = {r["id"]: r for r in old["cases"]} if old else {}
    recs = []
    for c in cases:
        rr = results.get(c["id"], {})
        base = prev.get(c["id"], {})
        fr = rr.get("fricas", base.get("fricas", {"status": "SKIPPED"}))
        rec = dict(c, fricas=fr, methods=dict(base.get("methods", {})))
        for mt in methods:
            rec["methods"][mt] = rr[mt]
        for mt, m in rec["methods"].items():      # FriCAS may have changed: reclassify all
            m["class"] = classify(m, fr, c["expect"])
        recs.append(rec)
    if old:
        done_ids = {r["id"] for r in recs}
        recs = [next(r for r in recs if r["id"] == o["id"]) if o["id"] in done_ids else o
                for o in old["cases"]]
        meta = dict(old, generated=time.strftime("%Y-%m-%d %H:%M:%S"),
                    merged=old.get("merged", []) + [{"methods": methods, "cases": len(cases),
                                                     "fricas": bool(fbin),
                                                     "elapsed": round(time.time() - t0, 1)}])
        meta["cases"] = recs
    else:
        meta = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
                "elapsed": round(time.time() - t0, 1),
                "timeout": args.timeout, "fricas_timeout": args.fricas_timeout,
                "methods": methods, "params": PARAMS, "cases": recs}
    with open(args.out, "w") as fh:
        json.dump(meta, fh, indent=1)
    methods = meta["methods"]
    for mt in methods:
        cnt = {}
        for r in recs:
            cnt[r["methods"][mt]["class"]] = cnt.get(r["methods"][mt]["class"], 0) + 1
        print(mt, cnt, file=sys.stderr)
    print(f"wrote {args.out} in {time.time() - t0:.0f} s", file=sys.stderr)


if __name__ == "__main__":
    main()
