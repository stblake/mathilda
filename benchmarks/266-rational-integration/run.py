#!/usr/bin/env python3.11
"""run.py -- runs benchmark 266 (rational integration) in Mathilda, the local
reference system (wolframscript) and sympy (ratint), and writes one TSV of raw
rows per system under results/<label>/.

    python3.11 benchmarks/266-rational-integration/run.py --label before
    python3.11 benchmarks/266-rational-integration/run.py --label after --systems mathilda
    python3.11 benchmarks/266-rational-integration/run.py --label x --families rootsum,stress -j 6
    python3.11 benchmarks/266-rational-integration/run.py --label grid --systems grid --timeout 15 -j 6

Resumable: a case already present in results/<label>/<system>.tsv is skipped,
so an interrupted run continues where it stopped (pass --fresh to start over).

Mathilda: one process per (case, method) -- a hang or a crash costs one case,
never the batch.  Each case is TimeConstrained (--timeout) inside the process and
the process itself is killed at 2x that wall time or at --rss-gb resident memory
(TimeConstrained cannot interrupt a single GMP/FLINT call).  The answer is
checked by PossibleZeroQ[D[F, x] - f] and, for a real-coefficient integrand free
of parameters, by a numeric continuity scan of Re F on [-4, 4] (a jump that the
integrand's poles do not explain is a spurious discontinuity, Bronstein 2.8).

Reference: one wolframscript process for the whole corpus; time and LeafCount
only.  It is never a correctness oracle.

sympy: ratint(f, x) in a child process per case with the same timeout; LeafCount
computed with the reference's counting rules (Rational = 3 leaves, I = 3, ...).
"""
import argparse
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CORPUS = os.path.join(HERE, "corpus.tsv")
WOLFRAM = os.environ.get("WOLFRAMSCRIPT", "wolframscript")
PY = os.environ.get("SYMPY_PYTHON", "python3.11")


def load_corpus(families=None, ids=None):
    rows = []
    with open(CORPUS) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cid, fam, wl, sy, note = (line.rstrip("\n").split("\t") + [""] * 5)[:5]
            if families and fam not in families:
                continue
            if ids and cid not in ids:
                continue
            rows.append(dict(id=cid, family=fam, wl=wl, sympy=sy, note=note))
    return rows


def done_ids(path):
    out = set()
    if os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                p = line.split("\t")
                if len(p) > 2:
                    out.add((p[0], p[1]))
    return out


# ------------------------------------------------------------------ Mathilda
M_DRIVER = r'''
f = %(f)s;
$tl = %(tl)s;
r = TimeConstrained[MemoryConstrained[AbsoluteTiming[%(call)s], %(mem)s, $MemOut], $tl, $TimedOut];
row = Which[
  r === $TimedOut, {"TIMEOUT", $tl, "-", "-", "-"},
  r === $MemOut, {"MEMOUT", "-", "-", "-", "-"},
  !FreeQ[r[[2]], Integrate], {"UNEVAL", r[[1]], LeafCount[r[[2]]], "-", "-"},
  True, Module[{F = r[[2]], v, jump = "-", real},
    v = TimeConstrained[PossibleZeroQ[D[F, x] - f], 60, $TimedOut];
    real = FreeQ[f, Complex | Root] && Complement[Variables[f], {x}] === {} && FreeQ[f, Power[_, _Rational]];
    If[real,
      Module[{h = 1/25, pts, Fv, fv, bad = 0, poles},
        pts = Range[-4, 4, h];
        Fv = Quiet[Re[N[F /. x -> #]] & /@ pts];
        fv = Quiet[N[f /. x -> #] & /@ (pts + h/2)];
        (* the real poles of f: an interval holding one is skipped (a log
           singularity of small residue is not a spurious jump) *)
        poles = Quiet[Re[Select[N[x /. NSolve[Denominator[Together[f]] == 0, x]], Abs[Im[#]] < 10^-6 &]]];
        If[!ListQ[poles], poles = {}];
        Do[
          If[NumberQ[Fv[[i]]] && NumberQ[Fv[[i + 1]]] && NumberQ[fv[[i]]] &&
             (poles === {} || Min[Abs[poles - N[pts[[i]] + h/2]]] > h) &&
             Abs[fv[[i]]] < 20 && Abs[N[f /. x -> pts[[i]]]] < 20 && Abs[N[f /. x -> pts[[i + 1]]]] < 20 &&
             Abs[Fv[[i + 1]] - Fv[[i]] - h fv[[i]]] > 3/10, bad++],
          {i, Length[pts] - 1}];
        jump = If[bad > 0, "JUMP" <> ToString[bad], "cont"]]];
    {Which[v === True, "PASS", v === $TimedOut, "VTIMEOUT", True, "FAIL"], r[[1]], LeafCount[F], jump,
     ToString[InputForm[F]]}]];
Print["ROW\t", StringRiffle[ToString /@ row, "\t"]];
'''


def mathilda_rss_kb(pid):
    try:
        out = subprocess.run(["ps", "-o", "rss=", "-p", str(pid)], capture_output=True, text=True).stdout
        return int(out.strip() or 0)
    except Exception:
        return 0


def run_mathilda_case(case, method, args):
    call = "Integrate[f, x]" if method == "Automatic" else f'Integrate[f, x, Method -> "{method}"]'
    src = M_DRIVER % dict(f=case["wl"], tl=args.timeout, call=call, mem=int(args.rss_gb * 0.8 * 2**30))
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False, dir=args.tmp) as fh:
        fh.write(src)
        path = fh.name
    t0 = time.time()
    p = subprocess.Popen([args.mathilda, "-file", path], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         text=True, cwd=ROOT)
    killed = None
    wall = 2 * args.timeout + 75
    while p.poll() is None:
        time.sleep(0.2)
        if time.time() - t0 > wall:
            killed = "KILLED_WALL"
        elif mathilda_rss_kb(p.pid) > args.rss_gb * 2**20:
            killed = "KILLED_RSS"
        if killed:
            p.kill()
            break
    out, err = p.communicate()
    os.unlink(path)
    row = None
    for line in out.splitlines():
        if line.startswith("ROW\t"):
            row = line[4:].split("\t", 4)
    if row is None:
        status = killed or ("CRASH" if p.returncode else "NOROW")
        row = [status, "-", "-", "-", (err.strip().splitlines() or ["-"])[-1][:200]]
    msgs = sum(1 for l in err.splitlines() if "::" in l)
    return [case["id"], method] + row[:4] + [str(msgs), row[4] if len(row) > 4 else "-"]


def run_mathilda(cases, args, outdir):
    path = os.path.join(outdir, "mathilda.tsv")
    have = done_ids(path)
    todo = [(c, m) for c in cases for m in args.methods if (c["id"], m) not in have]
    print(f"mathilda: {len(todo)} runs ({len(have)} cached)", file=sys.stderr)
    with open(path, "a") as fh, ThreadPoolExecutor(args.jobs) as ex:
        n = 0
        for row in ex.map(lambda cm: run_mathilda_case(cm[0], cm[1], args), todo):
            fh.write("\t".join(str(v).replace("\t", " ").replace("\n", " ") for v in row) + "\n")
            fh.flush()
            n += 1
            if n % 20 == 0 or row[2] not in ("PASS",):
                print(f"  [{n}/{len(todo)}] {row[0]} {row[1]} {row[2]} {row[3]}s lc={row[4]}", file=sys.stderr)


# ------------------------------------------------------------------ reference
W_DRIVER = r'''
cases = Get["%(cases)s"];
Do[
  Module[{id = bench266case[[1]], f = bench266case[[2]], r},
    r = TimeConstrained[AbsoluteTiming[Integrate[f, x]], %(tl)s, $TimedOut];
    row = Which[r === $TimedOut, {"TIMEOUT", %(tl)s, "-"},
       !FreeQ[r[[2]], Integrate], {"UNEVAL", r[[1]], LeafCount[r[[2]]]},
       True, {"OK", r[[1]], LeafCount[r[[2]]]}];
    WriteString["%(out)s", StringRiffle[Join[{id, "Automatic"}, ToString /@ row,
       {ToString[InputForm[If[r === $TimedOut, "-", r[[2]]]], PageWidth -> Infinity]}], "\t"] <> "\n"]],
  {bench266case, cases}];
Close["%(out)s"];
'''


def run_wolfram(cases, args, outdir):
    path = os.path.join(outdir, "wolfram.tsv")
    have = done_ids(path)
    todo = [c for c in cases if (c["id"], "Automatic") not in have]
    print(f"wolfram: {len(todo)} runs ({len(have)} cached)", file=sys.stderr)
    if not todo:
        return
    if not os.path.exists(WOLFRAM):
        print("  wolframscript not found; skipped", file=sys.stderr)
        return
    cpath = os.path.join(args.tmp, "w_cases.m")
    with open(cpath, "w") as fh:
        fh.write("{" + ",\n".join(f'{{"{c["id"]}", {c["wl"]}}}' for c in todo) + "}\n")
    part = os.path.join(args.tmp, "w_out.tsv")
    if os.path.exists(part):
        os.unlink(part)
    dpath = os.path.join(args.tmp, "w_driver.wl")
    with open(dpath, "w") as fh:
        fh.write(W_DRIVER % dict(cases=cpath, tl=args.timeout, out=part))
    subprocess.run([WOLFRAM, "-file", dpath], cwd=HERE)
    if os.path.exists(part):
        with open(path, "a") as fh, open(part) as ph:
            fh.write(ph.read())


# ------------------------------------------------------------------ sympy
SYMPY_CHILD = r'''
import sys, time
from sympy import *
from sympy.integrals.rationaltools import ratint
x, a, b, c = symbols("x a b c")

def lc(e):
    if e.is_Integer or e.is_Float or e.is_Symbol:
        return 1
    if e.is_Rational:
        return 3
    if e is S.ImaginaryUnit:
        return 3
    if isinstance(e, (Lambda,)):
        return 1 + sum(lc(t) for t in e.args[1:]) + len(e.args[0])
    try:
        args = e.args
    except Exception:
        return 1
    if not args:
        return 1
    return 1 + sum(lc(t) for t in args if isinstance(t, Basic))

f = sympify(sys.argv[1], locals=dict(x=x, a=a, b=b, c=c))
t0 = time.perf_counter()
F = ratint(f, x)
t = time.perf_counter() - t0
print("ROW\tOK\t%.6f\t%d\t%s" % (t, lc(F), str(F).replace("\t", " ").replace("\n", " ")))
'''


def run_sympy_case(case, args):
    if case["sympy"] == "-":
        return [case["id"], "ratint", "NA", "-", "-", "untranslatable (Root object)"]
    try:
        p = subprocess.run([PY, "-c", SYMPY_CHILD, case["sympy"]], capture_output=True, text=True,
                           timeout=args.timeout + 10)
    except subprocess.TimeoutExpired:
        return [case["id"], "ratint", "TIMEOUT", str(args.timeout), "-", "-"]
    for line in p.stdout.splitlines():
        if line.startswith("ROW\t"):
            _, st, t, l, F = line.split("\t", 4)
            return [case["id"], "ratint", st, t, l, F[:20000]]
    last = (p.stderr.strip().splitlines() or ["?"])[-1]
    return [case["id"], "ratint", "ERROR", "-", "-", last[:200]]


def run_sympy(cases, args, outdir):
    path = os.path.join(outdir, "sympy.tsv")
    have = done_ids(path)
    todo = [c for c in cases if (c["id"], "ratint") not in have]
    print(f"sympy: {len(todo)} runs ({len(have)} cached)", file=sys.stderr)
    with open(path, "a") as fh, ThreadPoolExecutor(args.jobs) as ex:
        for row in ex.map(lambda c: run_sympy_case(c, args), todo):
            fh.write("\t".join(row) + "\n")
            fh.flush()


# ------------------------------------------------------------------ algorithm grid
# Every combination of the rational integrator's algorithm choices
# (Method -> {"BronsteinRational", "Hermite" -> h, "LogPart" -> l, "RealForm" -> r}),
# all of them in one Mathilda process per case; each distinct answer is verified
# (PossibleZeroQ and the continuity scan) once.
GRID_HERMITE = ["Linear", "Quadratic", "HorowitzOstrogradsky"]
GRID_LOGPART = ["LazardRiobooTrager", "RothsteinTrager", "Czichowski"]
GRID_REALFORM = ["Automatic", "\"Rioboo\"", "\"LogToReal\"", "None"]

G_DRIVER = r'''
f = %(f)s;
$tl = %(tl)s;
seen266 = {};
check266[F_] := Module[{hit = Select[seen266, #[[1]] === F &, 1], v, jump = "-", real},
  If[hit =!= {}, Return[hit[[1, 2]]]];
  v = TimeConstrained[PossibleZeroQ[D[F, x] - f], 60, $TimedOut];
  real = FreeQ[f, Complex | Root] && Complement[Variables[f], {x}] === {} && FreeQ[f, Power[_, _Rational]];
  If[real,
    Module[{h = 1/25, pts, Fv, fv, bad = 0, poles},
      pts = Range[-4, 4, h];
      Fv = Quiet[Re[N[F /. x -> #]] & /@ pts];
      fv = Quiet[N[f /. x -> #] & /@ (pts + h/2)];
      poles = Quiet[Re[Select[N[x /. NSolve[Denominator[Together[f]] == 0, x]], Abs[Im[#]] < 10^-6 &]]];
      If[!ListQ[poles], poles = {}];
      Do[If[NumberQ[Fv[[i]]] && NumberQ[Fv[[i + 1]]] && NumberQ[fv[[i]]] &&
            (poles === {} || Min[Abs[poles - N[pts[[i]] + h/2]]] > h) &&
            Abs[fv[[i]]] < 20 && Abs[N[f /. x -> pts[[i]]]] < 20 && Abs[N[f /. x -> pts[[i + 1]]]] < 20 &&
            Abs[Fv[[i + 1]] - Fv[[i]] - h fv[[i]]] > 3/10, bad++], {i, Length[pts] - 1}];
      jump = If[bad > 0, "JUMP" <> ToString[bad], "cont"]]];
  hit = {Which[v === True, "PASS", v === $TimedOut, "VTIMEOUT", True, "FAIL"], jump};
  AppendTo[seen266, {F, hit}];
  hit];
Do[
  Module[{r = TimeConstrained[AbsoluteTiming[Integrate[f, x, Method -> Join[{"BronsteinRational"}, cfg[[2]]]]], $tl, $TimedOut], row},
    row = Which[
      r === $TimedOut, {"TIMEOUT", $tl, "-", "-"},
      !FreeQ[r[[2]], Integrate], {"UNEVAL", r[[1]], LeafCount[r[[2]]], "-"},
      True, Module[{c = check266[r[[2]]]}, {c[[1]], r[[1]], LeafCount[r[[2]]], c[[2]]}]];
    Print["GRID\t", cfg[[1]], "\t", StringRiffle[ToString /@ row, "\t"]]],
  {cfg, %(cfgs)s}];
'''


def grid_configs():
    out = []
    for h in GRID_HERMITE:
        for l in GRID_LOGPART:
            for r in GRID_REALFORM:
                name = f"{h}|{l}|{r.strip(chr(34))}"
                rule = f'{{"{name}", {{"Hermite" -> "{h}", "LogPart" -> "{l}", "RealForm" -> {r}}}}}'
                out.append((name, rule))
    return out


def run_grid_case(case, args, cfgs):
    src = G_DRIVER % dict(f=case["wl"], tl=args.timeout, cfgs="{" + ", ".join(r for _, r in cfgs) + "}")
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False, dir=args.tmp) as fh:
        fh.write(src)
        path = fh.name
    t0 = time.time()
    p = subprocess.Popen([args.mathilda, "-file", path], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                         text=True, cwd=ROOT)
    killed = None
    wall = len(cfgs) * (args.timeout + 5) + 120
    while p.poll() is None:
        time.sleep(0.2)
        if time.time() - t0 > wall:
            killed = "KILLED_WALL"
        elif mathilda_rss_kb(p.pid) > args.rss_gb * 2**20:
            killed = "KILLED_RSS"
        if killed:
            p.kill()
            break
    out, _ = p.communicate()
    os.unlink(path)
    rows, done = [], set()
    for line in out.splitlines():
        if line.startswith("GRID\t"):
            parts = line.split("\t")
            rows.append([case["id"]] + parts[1:6])
            done.add(parts[1])
    for name, _ in cfgs:
        if name not in done:
            rows.append([case["id"], name, killed or "CRASH", "-", "-", "-"])
    return rows


def run_grid(cases, args, outdir):
    path = os.path.join(outdir, "grid.tsv")
    have = {i for (i, _) in done_ids(path)}
    todo = [c for c in cases if c["id"] not in have]
    cfgs = grid_configs()
    print(f"grid: {len(todo)} cases x {len(cfgs)} configurations ({len(have)} cached)", file=sys.stderr)
    with open(path, "a") as fh, ThreadPoolExecutor(args.jobs) as ex:
        n = 0
        for rows in ex.map(lambda c: run_grid_case(c, args, cfgs), todo):
            for row in rows:
                fh.write("\t".join(str(v) for v in row) + "\n")
            fh.flush()
            n += 1
            bad = [r for r in rows if r[2] != "PASS" or r[5].startswith("JUMP")]
            if n % 25 == 0 or bad:
                print(f"  [{n}/{len(todo)}] {rows[0][0]}: {len(bad)} non-PASS/JUMP", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--label", required=True, help="results/<label>/ (e.g. before, after)")
    ap.add_argument("--systems", default="mathilda,wolfram,sympy")
    ap.add_argument("--methods", default="Automatic,BronsteinRational")
    ap.add_argument("--families", default="")
    ap.add_argument("--ids", default="")
    ap.add_argument("--timeout", type=int, default=30, help="seconds per integral")
    ap.add_argument("--rss-gb", type=float, default=4.0)
    ap.add_argument("-j", "--jobs", type=int, default=4)
    ap.add_argument("--mathilda", default=os.path.join(ROOT, "Mathilda"))
    ap.add_argument("--fresh", action="store_true")
    args = ap.parse_args()
    args.methods = args.methods.split(",")
    fams = set(filter(None, args.families.split(",")))
    ids = set(filter(None, args.ids.split(",")))
    cases = load_corpus(fams, ids)
    outdir = os.path.join(HERE, "results", args.label)
    os.makedirs(outdir, exist_ok=True)
    args.tmp = tempfile.mkdtemp(prefix="b266_")
    systems = args.systems.split(",")
    if args.fresh:
        for s in systems:
            p = os.path.join(outdir, f"{s}.tsv")
            if os.path.exists(p):
                os.unlink(p)
    t0 = time.time()
    # the reference and sympy run alongside Mathilda; each is single-threaded
    # (sympy) or one process (reference), so they do not crowd the Mathilda jobs
    with ThreadPoolExecutor(3) as ex:
        futs = []
        if "wolfram" in systems:
            futs.append(ex.submit(run_wolfram, cases, args, outdir))
        if "sympy" in systems:
            futs.append(ex.submit(run_sympy, cases, args, outdir))
        if "mathilda" in systems:
            futs.append(ex.submit(run_mathilda, cases, args, outdir))
        if "grid" in systems:
            futs.append(ex.submit(run_grid, cases, args, outdir))
        for f in futs:
            f.result()
    with open(os.path.join(outdir, "meta.txt"), "a") as fh:
        ver = subprocess.run([args.mathilda, "--version"], capture_output=True, text=True).stdout.strip()
        fh.write(f"{time.strftime('%Y-%m-%d %H:%M')} systems={args.systems} methods={','.join(args.methods)} "
                 f"timeout={args.timeout} jobs={args.jobs} cases={len(cases)} wall={time.time() - t0:.0f}s\n"
                 f"  {ver.splitlines()[0] if ver else '?'}\n")
    print(f"done in {time.time() - t0:.0f}s -> {outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
