#!/usr/bin/env python3
"""param_ab.py -- the Phase 4.4 step-6 A/B run (benchmark 258): does the
Automatic Integrate cascade get any worse when its Trager and Bronstein stages
accept integrands with symbolic parameters?

Corpus: the algebraic (A), mixed (M) and radical-at-the-bottom (B) cases of
benchmarks 256 and 257, each with one integer coefficient k turned into k + a
(param_gen.m): an `elem` case is parametrised in its antiderivative, so it stays
elementary over Q(a); the others are parametrised in the integrand (`unknown`).

Runs, per case: A = Integrate[f, x] with Integrate`$CascadeParameters = False
(no parametric input reaches either stage, the state before step 6), B = the
same with the value --b-flag (default True: both stages; the shipped default
{"Trager"} is the run that decided it, see benchmarks/258-parametric-ab/README.md).  Each
run verifies its own answer at the case's three points with a = 1/7 (40 digits,
|D[F] - f| <= 10^-20 max(1, |f|)), and records Integrate::nonelem as a
certificate -- the shootout's protocol (adapters.py), with the parameter value.

Gate (BRONSTEIN_PLAN.md, Phase 4.4 step 6): no WRONG, no FALSE_CERT (a CERT on
an `elem` case), no DISAGREE (B certifies what A answers), and no case where B
is worse than A: B loses A's verified answer, or its ELEM answer is much larger
(> 1.5 LeafCount + 20) or much slower (> 2 t + 5 s).

    python3 tools/integration_shootout/param_ab.py [--jobs 6] [--timeout 30] [--limit N]
        [--b-flag '{"Trager"}'] [--stem results-trager] [--rerun]
"""
import argparse
import json
import os
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(__file__))
import adapters as AD   # noqa: E402
import shootout as SH   # noqa: E402

REPO = AD.REPO
BENCH = os.path.join(REPO, "benchmarks", "258-parametric-ab")
SOURCES = [os.path.join(REPO, "benchmarks", b, "corpus.tsv")
           for b in ("256-integration-shootout", "257-compact-integrals")]
A0 = "1/7"

AB_TMPL = r"""
abRun[id_String, f_, pts_List, tmo_, flag_] := Module[{t0, res, cert = False, st, nums = {}, txt, t, r0, garbage},
  t0 = AbsoluteTime[];
  res = Block[{Integrate`$CascadeParameters = flag},
    TimeConstrained[Check[r0 = Integrate[f, x], cert = True; r0, {Integrate::nonelem}], tmo, $Aborted]];
  t = AbsoluteTime[] - t0;
  st = Which[res === $Aborted, "TIMEOUT", !FreeQ[res, Integrate], "UNEVAL", True, "ANSWER"];
  garbage = st === "ANSWER" && (! FreeQ[res, $Failed | $Aborted] ||
    ! FreeQ[res, s_Symbol /; StringContainsQ[Context[s], "Private`"]]);
  If[garbage, nums = {"bad"}];
  If[st === "ANSWER" && ! garbage,
    nums = Table[Module[{d, fv},
        d = TimeConstrained[Quiet[N[((D[res, x] - f) /. a -> A0) /. x -> p, 40]], VERIFY, $Aborted];
        fv = Quiet[N[(f /. a -> A0) /. x -> p, 40]];
        Which[!NumberQ[d] || !NumberQ[fv], "na",
              Abs[d] <= 10^-TOL Max[1, Abs[fv]], "ok", True, "bad"]], {p, pts}]];
  txt = If[res === $Aborted, "", ToString[res, InputForm]];
  txt = StringReplace[txt, {"\t" -> " ", "\n" -> " "}];
  Print["@@R\t", id, "\t", st, "\t", ToString[N[t]], "\t", StringRiffle[nums, ","], "\t",
    ToString[cert], "\t", ToString[LeafCount[res]] <> "|" <> StringTake[txt, UpTo[TRUNC]]]];
"""


def gen_corpus(limit=None):
    """The parametric corpus (param_gen.m over the A/M/B cases of 256 and 257)."""
    src = []
    for path in SOURCES:
        for c in SH.load_corpus(path):
            if c["class"] in ("A", "M", "B"):
                src.append(c)
    if limit:
        src = src[::max(1, len(src) // limit)][:limit]
    fd, data = tempfile.mkstemp(suffix=".m")
    with os.fdopen(fd, "w") as fh:
        fh.write("pgCases = {\n")
        fh.write(",\n".join('{"%s", "%s", Hold[%s], %s}' % (
            c["id"], c["expect"], c["f"], "None" if c["F"] is None else "Hold[%s]" % c["F"]) for c in src))
        fh.write("};\n")
        fh.write('Get["%s"];\n' % os.path.join(REPO, "tools", "integration_shootout", "param_gen.m"))
    out, err, to, wall = AD.run_proc([AD.MATHILDA, "-file", data], 3600 + 30 * len(src))
    os.unlink(data)
    cases = []
    for ln in out.splitlines():
        p = ln.split("\t")
        if len(p) == 5 and p[0] == "PCASE":
            cases.append({"id": p[1], "expect": p[2], "f": p[3], "points": p[4]})
    return cases


def run_batch(cases, flag, timeout):
    head = (AB_TMPL.replace("A0", A0).replace("VERIFY", str(AD.VERIFY))
            .replace("TOL", str(AD.TOL_EXP)).replace("TRUNC", str(AD.TRUNC)))
    body = "".join(f'abRun["{c["id"]}", {c["f"]}, {c["points"]}, {timeout}, {flag}];\n' for c in cases)
    fd, path = tempfile.mkstemp(suffix=".m")
    with os.fdopen(fd, "w") as fh:
        fh.write(head + body)
    out, err, to, wall = AD.run_proc([AD.MATHILDA, "-file", path], 20 + len(cases) * (timeout + 3 * AD.VERIFY + 5))
    os.unlink(path)
    recs = AD.parse_rlines(out)
    for r in recs.values():
        lc, _, txt = r["text"].partition("|")
        r["leaf"] = int(lc) if lc.isdigit() else None
        r["text"] = txt
    return recs


def run_pairs(cases, timeout, jobs, bflag="True", batch=6):
    """A and B per batch, back to back in one worker: both arms of a case see the
    same machine load, so a timing difference is the stage's, not the neighbours'."""
    chunks = [cases[i:i + batch] for i in range(0, len(cases), batch)]
    ra, rb = {}, {}

    def arm(ch, flag):
        got = run_batch(ch, flag, timeout)
        for c in ch:            # a case missing from its batch: alone
            if c["id"] not in got:
                got.update(run_batch([c], flag, timeout))
            if c["id"] not in got:
                got[c["id"]] = {"id": c["id"], "status": "ERROR", "time": None, "nums": [], "cert": False,
                                "text": "", "leaf": None}
        return got
    with ThreadPoolExecutor(jobs) as ex:
        for ga, gb in ex.map(lambda ch: (arm(ch, "False"), arm(ch, bflag)), chunks):
            ra.update(ga)
            rb.update(gb)
    return ra, rb


def verdict(raw):
    return SH.classify_run("mathilda", raw)


def compare(c, ra, rb):
    """The gate's findings for one case (hard, worse, better)."""
    va, vb = verdict(ra), verdict(rb)
    hard, worse, better = [], [], []
    if vb == "WRONG":
        hard.append("WRONG")
    if vb == "CERT" and c["expect"] == "elem":
        hard.append("FALSE_CERT")
    if vb == "CERT" and va in ("ELEM", "SPECIAL"):
        hard.append("DISAGREE")
    if va in ("ELEM", "SPECIAL") and vb not in ("ELEM", "SPECIAL"):
        worse.append(f"lost {va} -> {vb}")
    if va == "ELEM" and vb == "SPECIAL":
        worse.append("ELEM -> SPECIAL")
    if va == vb == "ELEM" and ra.get("leaf") and rb.get("leaf") and rb["leaf"] > 1.5 * ra["leaf"] + 20:
        worse.append(f"larger {ra['leaf']} -> {rb['leaf']}")
    if va == vb and ra.get("time") is not None and rb.get("time") is not None and rb["time"] > 2 * ra["time"] + 5:
        worse.append(f"slower {ra['time']:.1f}s -> {rb['time']:.1f}s")
    if va not in ("ELEM", "SPECIAL") and vb in ("ELEM", "SPECIAL"):
        better.append(f"{va} -> {vb}")
    if va in ("UNEVAL", "TIMEOUT") and vb == "CERT":
        better.append(f"{va} -> CERT")
    if va == vb == "ELEM" and ra.get("leaf") and rb.get("leaf") and rb["leaf"] * 1.5 + 20 < ra["leaf"]:
        better.append(f"smaller {ra['leaf']} -> {rb['leaf']}")
    return va, vb, hard, worse, better


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--timeout", type=int, default=30)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--fresh", action="store_true", help="regenerate the parametric corpus")
    ap.add_argument("--b-flag", default="True",
                    help='Integrate`$CascadeParameters in arm B: True, {"Trager"}, {"Bronstein"}')
    ap.add_argument("--stem", default=None,
                    help="output name: <stem>.jsonl and <STEM>.md (default results/REPORT, rerun/RERUN)")
    ap.add_argument("--rerun", action="store_true",
                    help="re-run only the cases results.jsonl flags (hard or worse), A and B again, "
                         "and write rerun.jsonl/RERUN.md: a timing or timeout finding made under load "
                         "is confirmed or dropped on a quiet machine")
    args = ap.parse_args()
    os.makedirs(BENCH, exist_ok=True)
    cpath = os.path.join(BENCH, "corpus.jsonl")
    if args.fresh or not os.path.exists(cpath):
        cases = gen_corpus(args.limit)
        with open(cpath, "w") as fh:
            for c in cases:
                fh.write(json.dumps(c) + "\n")
    cases = [json.loads(ln) for ln in open(cpath)]
    if args.rerun:
        flagged = {r["id"] for r in map(json.loads, open(os.path.join(BENCH, (args.stem or "results") + ".jsonl")))
                   if r["hard"] or r["worse"]}
        cases = [c for c in cases if c["id"] in flagged]
    if args.limit:
        cases = cases[:args.limit]
    print(f"{len(cases)} parametric cases", flush=True)
    ra, rb = run_pairs(cases, args.timeout, args.jobs, args.b_flag)
    rows, counts = [], {}
    for c in cases:
        va, vb, hard, worse, better = compare(c, ra[c["id"]], rb[c["id"]])
        counts[(va, vb)] = counts.get((va, vb), 0) + 1
        rows.append({"id": c["id"], "expect": c["expect"], "f": c["f"], "A": va, "B": vb,
                     "tA": ra[c["id"]].get("time"), "tB": rb[c["id"]].get("time"),
                     "leafA": ra[c["id"]].get("leaf"), "leafB": rb[c["id"]].get("leaf"),
                     "hard": hard, "worse": worse, "better": better})
    stem = ("rerun-" + args.stem if args.stem else "rerun") if args.rerun else (args.stem or "results")
    with open(os.path.join(BENCH, stem + ".jsonl"), "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    hard = [r for r in rows if r["hard"]]
    worse = [r for r in rows if r["worse"]]
    better = [r for r in rows if r["better"]]
    lines = ["# 258 — parametric A/B of the Trager and Bronstein cascade stages", "",
             f"{len(cases)} cases; A = `Integrate`$CascadeParameters = False`, "
             f"B = `{args.b_flag}`; timeout {args.timeout} s, check at a = {A0}.", "",
             f"- hard failures in B: **{len(hard)}**", f"- B worse than A: **{len(worse)}**",
             f"- B better than A: **{len(better)}**", "", "| A | B | cases |", "|---|---|---|"]
    for (va, vb), n in sorted(counts.items(), key=lambda kv: -kv[1]):
        lines.append(f"| {va} | {vb} | {n} |")
    for title, sel, key in (("Hard failures", hard, "hard"), ("Worse", worse, "worse"), ("Better", better, "better")):
        lines += ["", f"## {title}", ""]
        for r in sel:
            lines.append(f"- `{r['id']}` ({r['expect']}): {', '.join(r[key])} — `{r['f'][:160]}`")
    with open(os.path.join(BENCH, ("REPORT" if stem == "results" else stem.upper()) + ".md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines[:14]))
    return 1 if hard or worse else 0


if __name__ == "__main__":
    sys.exit(main())
