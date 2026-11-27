#!/usr/bin/env python3
"""goursat_stress.py -- benchmark 265-goursat-stress: the GeneralisedGoursat corpora
through two Mathilda Integrate methods, goursat.py and plain sympy.

    python3 tools/goursat/goursat_stress.py                         # the timed run, 30 s per case
    python3 tools/goursat/goursat_stress.py --correctness           # wrong answers only: 60 s, 4 jobs
    python3 tools/goursat/goursat_stress.py --with-sympy            # add plain sympy integrate (slow)
    python3 tools/goursat/goursat_stress.py --corpora stress_7 --methods GeneralisedGoursat
    python3 tools/goursat/goursat_stress.py --families sqrt-klein --fresh
    python3 tools/goursat/goursat_stress.py --goursatpy fresh       # re-time goursat.py (slow)
    python3 tools/goursat/goursat_stress.py --retabulate            # REPORT.md from the last run

CORPUS.  benchmarks/265-goursat-stress/corpus.tsv (tools/goursat/gen_stress_corpus.m):
one CASE row per case of the nine tests/data/generalised_goursat corpora, class G,
family = the corpus's Family, expect = its Truth, plus the corpus, label and truth.

RUNS (columns).
  Mathilda:Automatic           Integrate[f, x]
  Mathilda:GeneralisedGoursat  Integrate[f, x, Method -> "GeneralisedGoursat"] (the .m port),
                               with Integrate`$GeneralisedGoursatBudget = Infinity so the
                               per-case limit, not the stage's own 30 s budget, ends it.
                               A certified non-elementary integral comes back as a
                               decline or as partial + Inactive[Integrate][residual, x]
                               (UNEVAL either way), so the stage symbol
                               Integrate`GeneralisedGoursat[f, x] is asked again,
                               untimed, only to read the certificate.
  (The C GoursatAlgebraic column was dropped when that method was deleted, 2026-11-28;
  runs before then carry it in history.jsonl.)
  goursat.py                   the sympy reference.  By default its timings are READ from
                               tests/data/generalised_goursat/sympy_<corpus>.tsv
                               (tools/goursat/bench_sympy.py: goursat_integrate on F, R, p
                               or integrate_algebraic, reduction only, 120 s limit, its own
                               numerical verify).  --goursatpy fresh re-runs it through
                               adapters.run_goursatpy (integrate_algebraic on every case,
                               this script's timeout), and then its answers are judged by
                               PossibleZeroQ like everyone's.
  sympy                        plain sympy integrate (adapters.run_sympy); opt-in, --with-sympy:
                               the slowest column by far (half its runs time out at 30 s) and
                               the least informative.

MODES.  The default is the TIMED run (30 s, 6 jobs): run it on a quiet machine; its
REPORT.md is the baseline.  --correctness is the load-independent pass (60 s, 4 jobs by
default): the same runs and the same PossibleZeroQ verdicts, but REPORT.md is headed
"correctness pass; timings pending a quiet machine", its times are not a baseline, and
its timeouts are provisional.  Wrong answers and verdicts do not depend on load.

VERIFICATION.  Every Mathilda (and sympy) answer is judged by PossibleZeroQ[D[F, x] - f]
in Mathilda (crosscheck.run_batch, 60 s): zero -> verified, nonzero -> WRONG, anything
else -> unverified.  There is no numeric fallback for the verdict; each system's own
three-point check is recorded ("nums") and decides nothing.

PER-RUN CLASS: ELEM | SPECIAL (verified) | WRONG | UNVERIFIED | CERT (certified
non-elementary) | UNEVAL (declined / unevaluated) | TIMEOUT | ERROR.  Hard failures:
WRONG, and FALSE_CERT (a certificate on an expect-elem case).

CACHE.  results/cache.jsonl, keyed by system, method, case id, integrand and (Mathilda)
the binary + src/internal .m stamp, so an interrupted run resumes; --fresh ignores it.
The binary is copied to a temporary snapshot first (SHOOTOUT_MATHILDA pins another), so a
rebuild of ./Mathilda mid-run does not mix builds; the .m packages are read from
src/internal at run time and their stamp is recorded at the start and end of the run.
"""
import argparse
import collections
import concurrent.futures as cf
import datetime
import json
import os
import shutil
import statistics
import sys
import tempfile
import threading

os.environ.setdefault("SHOOTOUT_TRUNC", "1000000")   # full answers: PossibleZeroQ needs them
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "integration_shootout"))
import adapters as A  # noqa: E402
import crosscheck as X  # noqa: E402
import shootout as SH  # noqa: E402
import translate as T  # noqa: E402

BENCH = os.path.join(A.REPO, "benchmarks", "265-goursat-stress")
GGDATA = os.path.join(A.REPO, "tests", "data", "generalised_goursat")
METHODS = ["Automatic", "GeneralisedGoursat"]
# further Integrate methods that may be requested with --methods (not run by default)
EXTRA_METHODS = ["Trager", "TragerKauers", "Bronstein"]
CORPORA = ["unit", "examples", "literature", "welz", "stress_2026", "stress_7", "stress_11",
           "pseudoelliptic_sqrt", "pseudoelliptic_cbrt"]
PYLAB = "goursat.py"
CLOSED = ("ELEM", "SPECIAL")
DECIDED = ("ELEM", "SPECIAL", "CERT")
VERIFY = 10        # seconds per point of the own (recorded, non-deciding) check
MEM = 4 * 10 ** 9  # MemoryConstrained per case

# warm-up: every method's lazy loads (the .m packages, FLINT contexts) are paid once per
# process before the first timed case, so a batch's first case is not charged for them
WARM = {"Automatic": "Integrate[x/Sqrt[1 + x^4], x]; Integrate[1/(1 + x^3)^(1/3), x]",
        "GeneralisedGoursat": 'Integrate[(x^2 - 1)/((x^2 + 1) Sqrt[x^4 + 1]), x, Method -> "GeneralisedGoursat"]'}

GS_TMPL = r"""
Integrate`$GeneralisedGoursatBudget = Infinity;
gsWarm = TimeConstrained[Quiet[@WARM@], 60, $Aborted];
gsRun[id_String, f_, pts_List, tmo_] := Module[{t0, t, res, r0, r1, cert = False, st, nums = {}, txt,
    lc = "", ct = "", p0, garbage},
  t0 = AbsoluteTime[];
  res = MemoryConstrained[TimeConstrained[Check[r0 = @CALL@, cert = True; r0, {Integrate::nonelem}],
          tmo, $Aborted], @MEM@, "MEMORY"];
  t = AbsoluteTime[] - t0;
  st = Which[res === $Aborted, "TIMEOUT", res === "MEMORY", "MEMORY", !FreeQ[res, Integrate], "UNEVAL",
             True, "ANSWER"];
  (* the certificate the Method option drops: ask the stage itself, untimed *)
  If[@CERTPROBE@ && st === "UNEVAL" && !cert,
    p0 = AbsoluteTime[];
    r1 = MemoryConstrained[TimeConstrained[Quiet[Integrate`GeneralisedGoursat[f, x]], tmo, $Aborted],
           @MEM@, $Aborted];
    ct = ToString[N[AbsoluteTime[] - p0]];
    If[MatchQ[r1, {"not elementary", _, _}], cert = True]];
  (* $Failed / a package-private symbol in an "answer" is garbage reaching the user *)
  garbage = st === "ANSWER" && (! FreeQ[res, $Failed | $Aborted] ||
    ! FreeQ[res, s_Symbol /; StringContainsQ[Context[s], "Private`"]]);
  If[garbage, nums = {"bad"}];
  If[st === "ANSWER", lc = ToString[LeafCount[res]]];
  If[st === "ANSWER" && ! garbage,
    nums = Table[Module[{d, fv},
        d = TimeConstrained[Quiet[N[(D[res, x] - f) /. x -> p, 40]], @VERIFY@, $Aborted];
        fv = Quiet[N[f /. x -> p, 40]];
        Which[!NumberQ[d] || !NumberQ[fv], "na", Abs[d] <= 10^-20 Max[1, Abs[fv]], "ok", True, "bad"]],
      {p, pts}]];
  txt = If[st === "ANSWER" || st === "UNEVAL", ToString[res, InputForm], ""];
  txt = StringReplace[txt, {"\t" -> " ", "\n" -> " "}];
  Print["@@R\t", id, "\t", st, "\t", ToString[N[t]], "\t", StringRiffle[nums, ","], "\t",
    ToString[cert], "\t", lc, "\t", ct, "\t", StringTake[txt, UpTo[@TRUNC@]]]];
"""


def mathilda_call(method):
    return "Integrate[f, x]" if method == "Automatic" else f'Integrate[f, x, Method -> "{method}"]'


def parse_gs(out):
    """@@R id status time nums cert leafcount certtime text"""
    recs = {}
    for ln in out.splitlines():
        if not ln.startswith("@@R\t"):
            continue
        p = ln.split("\t", 8) + [""] * 9
        try:
            t = float(p[3])
        except ValueError:
            t = None
        st = p[2]
        detail = ""
        if st == "MEMORY":
            st, detail = "ERROR", "MemoryConstrained (4 GB)"
        recs[p[1]] = {"id": p[1], "status": st, "time": t, "nums": [s for s in p[4].split(",") if s],
                      "cert": p[5] == "True", "leaf": int(p[6]) if p[6].isdigit() else None,
                      "cert_time": float(p[7]) if p[7] else None, "text": p[8][:A.TRUNC],
                      "detail": detail}
    return recs


def run_mathilda(cases, method, timeout, tmpdir, batch=8):
    head = (GS_TMPL.replace("@WARM@", WARM.get(method, f'Integrate[(x^2 - 1)/((x^2 + 1) Sqrt[x^4 + 1]), x, Method -> "{method}"]')).replace("@CALL@", mathilda_call(method))
            .replace("@CERTPROBE@", "True" if method == "GeneralisedGoursat" else "False")
            .replace("@VERIFY@", str(VERIFY)).replace("@MEM@", str(MEM)).replace("@TRUNC@", str(A.TRUNC)))

    def runner(cs):
        body = "".join(f'gsRun["{c["id"]}", {c["f"]}, {c["points"]}, {timeout}];\n' for c in cs)
        fd, path = tempfile.mkstemp(suffix=".m", dir=tmpdir)
        with os.fdopen(fd, "w") as fh:
            fh.write(head + body)
        mult = 2 if method == "GeneralisedGoursat" else 1
        out, err, to, wall = A.run_proc([A.MATHILDA, "-file", path],
                                        60 + len(cs) * (mult * timeout + 3 * VERIFY + 10))
        os.unlink(path)
        return parse_gs(out), err
    return A.batched(runner, cases, timeout, batch)


# --------------------------------------------------------------------------- corpus
def load_cases(path):
    cases = SH.load_corpus(path)
    extra = {}
    for ln in open(path):
        p = ln.rstrip("\n").split("\t")
        if p[0] == "CASE" and len(p) >= 13:
            extra[p[1]] = {"corpus": p[10], "label": p[11], "truth": p[12]}
    for c in cases:
        c.update(extra.get(c["id"], {"corpus": "?", "label": c["id"], "truth": "None"}))
    return cases


PY_STATUS = {"elementary": "ELEM", "non-elementary": "CERT", "undecided": "UNEVAL",
             "TIMEOUT": "TIMEOUT", "ERROR": "ERROR"}


def goursatpy_from_tsv(cases):
    """goursat.py's verdict and time per case, from bench_sympy.py's sympy_<corpus>.tsv.
    Its own numerical verification (defect) is all there is: no antiderivative text."""
    by = {}
    for corp in CORPORA:
        path = os.path.join(GGDATA, f"sympy_{corp}.tsv")
        if not os.path.exists(path):
            continue
        for ln in list(open(path))[1:]:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 3:
                by[(corp, p[0])] = p
    out = {}
    for c in cases:
        p = by.get((c["corpus"], c["label"]))
        if p is None:
            continue
        cls = PY_STATUS.get(p[1], "ERROR")
        if cls == "ELEM" and p[3] and float(p[3]) > 1e-20:
            cls = "WRONG"
        out[c["id"]] = {"cls": cls, "time": float(p[2]), "source": "tsv", "own_defect": p[3],
                        "detail": p[5] if len(p) > 5 else ""}
    return out


# --------------------------------------------------------------------------- run
def classify(raw, xv):
    st = raw["status"]
    if st == "UNEVAL":
        return "CERT" if raw.get("cert") else "UNEVAL"
    if st != "ANSWER":
        return st
    if raw.get("nums") == ["bad"]:
        return "WRONG"                       # garbage: $Failed / Private` symbol
    if xv == "nonzero":
        return "WRONG"
    if xv != "zero":
        return "UNVERIFIED"
    return "SPECIAL" if T.special_q("mathilda", raw.get("text", "")) else "ELEM"


def internal_stamp():
    newest = 0
    for root, _d, files in os.walk(os.path.join(A.REPO, "src", "internal")):
        for fn in files:
            if fn.endswith(".m"):
                newest = max(newest, os.path.getmtime(os.path.join(root, fn)))
    return int(newest)


class Load:
    """os.getloadavg() sampled every 30 s while the run lasts"""
    def __init__(self):
        self.samples, self.stop = [os.getloadavg()[0]], threading.Event()
        self.th = threading.Thread(target=self.loop, daemon=True)
        self.th.start()

    def loop(self):
        while not self.stop.wait(30):
            self.samples.append(os.getloadavg()[0])

    def done(self):
        self.stop.set()
        s = self.samples + [os.getloadavg()[0]]
        return {"start": round(s[0], 1), "end": round(s[-1], 1), "mean": round(statistics.mean(s), 1),
                "max": round(max(s), 1), "cpus": os.cpu_count()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", default="mathilda,goursatpy")
    ap.add_argument("--with-sympy", action="store_true", help="add the plain sympy integrate column")
    ap.add_argument("--correctness", action="store_true",
                    help="load-independent correctness pass: 60 s per case, 4 jobs, report labelled as such")
    ap.add_argument("--methods", default=",".join(METHODS))
    ap.add_argument("--corpora", default=None, help="comma-separated corpus names")
    ap.add_argument("--families", default=None, help="comma-separated family prefixes")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--timeout", type=int, default=None, help="per case (default 30; 60 with --correctness)")
    ap.add_argument("--jobs", type=int, default=None, help="parallel processes (default 6; 4 with --correctness)")
    ap.add_argument("--goursatpy", choices=("tsv", "fresh"), default="tsv")
    ap.add_argument("--fresh", action="store_true")
    ap.add_argument("--bench", default=BENCH)
    ap.add_argument("--retabulate", action="store_true")
    ap.add_argument("--check", action="store_true", help="exit 1 on a WRONG or FALSE_CERT")
    a = ap.parse_args()
    if a.timeout is None:
        a.timeout = 60 if a.correctness else 30
    if a.jobs is None:
        a.jobs = 4 if a.correctness else 6
    if a.with_sympy and "sympy" not in a.systems.split(","):
        a.systems += ",sympy"
    bench = os.path.abspath(a.bench)
    rdir = os.path.join(bench, "results")
    os.makedirs(rdir, exist_ok=True)
    if a.retabulate:
        name = sorted(f[:-5] for f in os.listdir(rdir) if f[:4].isdigit() and f.endswith(".json"))[-1]
        saved = json.load(open(os.path.join(rdir, name + ".json")))
        write_report(os.path.join(bench, "REPORT.md"), saved, saved["date"])
        print(f"REPORT.md from results/{name}.json")
        return

    tmpdir = tempfile.mkdtemp(prefix="gstress-")
    if "SHOOTOUT_MATHILDA" not in os.environ:
        snap = os.path.join(tmpdir, "Mathilda")
        shutil.copy2(A.MATHILDA, snap)         # copy2 keeps the mtime, so the cache stamp is the build's
        A.MATHILDA = snap
    systems = [s for s in a.systems.split(",") if s]
    methods = [m for m in METHODS + EXTRA_METHODS if m in a.methods.split(",")]
    cases = load_cases(os.path.join(bench, "corpus.tsv"))
    if a.corpora:
        cases = [c for c in cases if c["corpus"] in a.corpora.split(",")]
    if a.families:
        cases = [c for c in cases if any(c["family"].startswith(f) for f in a.families.split(","))]
    if a.limit:
        cases = cases[:a.limit]
    cache = SH.Cache(os.path.join(rdir, "cache.jsonl"), a.fresh)
    stamp0, istamp0 = SH.mathilda_stamp(), internal_stamp()
    version = A.version("mathilda")
    load = Load()

    def ckey(s, m, c):
        return "|".join([s, m or "", c["id"], c["key"], stamp0 if s == "mathilda" else "", str(a.timeout)])

    units = []
    if "mathilda" in systems:
        for m in methods:
            units += [("mathilda", m, cases[i:i + 8]) for i in range(0, len(cases), 8)]
    if "sympy" in systems and A.available("sympy"):
        units += [("sympy", None, cases[i:i + 10]) for i in range(0, len(cases), 10)]
    if "goursatpy" in systems and a.goursatpy == "fresh" and A.goursatpy_available():
        units += [("goursatpy", None, cases[i:i + 10]) for i in range(0, len(cases), 10)]
    # longest-first: GeneralisedGoursat and sympy on the stress corpora dominate the wall
    units.sort(key=lambda u: (u[0] != "sympy", u[1] != "GeneralisedGoursat"))
    total, done, plock = sum(len(u[2]) for u in units), collections.Counter(), threading.Lock()

    def work(u):
        s, m, cs = u
        todo = [c for c in cs if cache.get(ckey(s, m, c)) is None]
        if todo:
            if s == "mathilda":
                got = run_mathilda(todo, m, a.timeout, tmpdir)
            elif s == "sympy":
                got = A.run_sympy(todo, a.timeout, tmpdir)
            else:
                got = A.run_goursatpy(todo, a.timeout, tmpdir)
            for c in todo:
                r = got.get(c["id"]) or {"id": c["id"], "status": "ERROR", "time": None, "nums": [],
                                          "cert": False, "text": "", "detail": "no record"}
                cache.put(ckey(s, m, c), r)
        with plock:
            done["n"] += len(cs)
            if done["n"] % 200 < len(cs):
                print(f"  {done['n']}/{total} runs", flush=True)

    t0 = datetime.datetime.now()
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        list(ex.map(work, units))

    # PossibleZeroQ cross-check of every answer
    runs_of = [("mathilda", m, f"Mathilda:{m}") for m in methods if "mathilda" in systems]
    if "sympy" in systems:
        runs_of.append(("sympy", None, "sympy"))
    if "goursatpy" in systems and a.goursatpy == "fresh":
        runs_of.append(("goursatpy", None, PYLAB))
    xv = crosscheck(cases, runs_of, cache, ckey, stamp0, tmpdir, a.jobs)
    wall = (datetime.datetime.now() - t0).total_seconds()
    loadinfo = load.done()

    pytsv = goursatpy_from_tsv(cases) if ("goursatpy" in systems and a.goursatpy == "tsv") else {}
    rows = []
    for c in cases:
        runs, raw = {}, {}
        for s, m, lab in runs_of:
            r = cache.get(ckey(s, m, c))
            if r is None:
                continue
            v = xv.get(ckey(s, m, c))
            runs[lab] = classify(r, v)
            raw[lab] = {k: r.get(k) for k in ("status", "time", "nums", "cert", "leaf", "cert_time", "detail")}
            if raw[lab]["leaf"] is None and r.get("status") == "ANSWER" and s != "mathilda":
                raw[lab]["leaf"] = None
            raw[lab]["pzq"] = v
            raw[lab]["text"] = (r.get("text") or "")[:400]
        if c["id"] in pytsv:
            p = pytsv[c["id"]]
            runs[PYLAB] = p["cls"]
            raw[PYLAB] = {"status": p["cls"], "time": p["time"], "source": "tsv",
                          "own_defect": p["own_defect"], "detail": p["detail"]}
        verd = []
        for lab, cls in runs.items():
            if lab.startswith("Mathilda:") and cls == "WRONG":
                verd.append(f"WRONG {lab}")
            if lab.startswith("Mathilda:") and cls == "CERT" and c["expect"] == "elem":
                verd.append(f"FALSE_CERT {lab}")
        rows.append({"id": c["id"], "corpus": c["corpus"], "family": c["family"], "expect": c["expect"],
                     "label": c["label"], "f": c["f"], "runs": runs, "raw": raw, "verdicts": verd})

    labels = [lab for _s, _m, lab in runs_of] + ([PYLAB] if pytsv else [])
    full = not (a.corpora or a.families or a.limit) and len(methods) == len(METHODS)
    date = datetime.date.today().isoformat()
    saved = {"date": date, "wall": wall, "timeout": a.timeout, "jobs": a.jobs, "load": loadinfo,
             "mathilda": version, "stamp": {"start": istamp0, "end": internal_stamp()},
             "goursatpy": a.goursatpy, "mode": "correctness" if a.correctness else "timed",
             "labels": labels, "rows": rows}
    if full:
        json.dump(saved, open(os.path.join(rdir, f"{date}{'-correctness' if a.correctness else ''}.json"), "w"),
                  indent=0)
        write_report(os.path.join(bench, "REPORT.md"), saved, date)
        with open(os.path.join(bench, "history.jsonl"), "a") as fh:
            fh.write(json.dumps(history_entry(saved)) + "\n")
        print(f"report: {os.path.join(bench, 'REPORT.md')}")
    else:
        # a partial run (a method or case subset) keeps its own JSON, for compare_pdf.py
        json.dump(saved, open(os.path.join(rdir, f"{date}-partial.json"), "w"), indent=0)
        write_report(os.path.join(bench, "REPORT-partial.md"), saved, date)
        print(f"report: {os.path.join(bench, 'REPORT-partial.md')}")
    hard = [(r["id"], v) for r in rows for v in r["verdicts"]]
    for cid, v in hard:
        print(f"*** {v}  [{cid}]")
    print(f"{len(rows)} cases, {wall / 60:.1f} min, load {loadinfo}")
    sys.exit(1 if (a.check and hard) else 0)


def crosscheck(cases, runs_of, cache, ckey, stamp, tmpdir, jobs):
    todo, out = [], {}
    for c in cases:
        for s, m, _lab in runs_of:
            k = ckey(s, m, c)
            r = cache.get(k)
            if r is None or r.get("status") != "ANSWER" or r.get("nums") == ["bad"]:
                continue
            xk = X.ckey(k + "|" + stamp, r.get("text") or "")
            hit = cache.get(xk)
            if hit is not None:
                out[k] = hit["verdict"]
                continue
            try:
                wl = r["text"] if s == "mathilda" else X.to_wl("sympy", r["text"])
            except X.Untranslatable as e:
                cache.put(xk, {"id": c["id"], "verdict": "untranslatable", "detail": str(e)})
                out[k] = "untranslatable"
                continue
            todo.append((k, xk, c, wl))
    if not todo:
        return out
    print(f"  cross-check: {len(todo)} answers to judge by PossibleZeroQ", flush=True)

    def work(chunk):
        got = X.run_batch([(f"x{j}", t[2]["fullform"], t[3]) for j, t in enumerate(chunk)], tmpdir)
        for j, (k, xk, c, wl) in enumerate(chunk):
            v = got.get(f"x{j}", "undecided")
            cache.put(xk, {"id": c["id"], "verdict": v, "detail": ""})
            out[k] = v
    with cf.ThreadPoolExecutor(max_workers=jobs) as ex:
        list(ex.map(work, [todo[i:i + 10] for i in range(0, len(todo), 10)]))
    return out


# --------------------------------------------------------------------------- report
def stats_of(rows, lab):
    cnt = collections.Counter(r["runs"][lab] for r in rows if lab in r["runs"])
    n = sum(cnt.values())
    dec_t = [r["raw"][lab]["time"] for r in rows if r["runs"].get(lab) in DECIDED and r["raw"][lab].get("time")]
    all_t = [r["raw"][lab]["time"] for r in rows if lab in r["raw"] and r["raw"][lab].get("time") is not None]
    leaf = [r["raw"][lab].get("leaf") for r in rows if r["runs"].get(lab) in CLOSED and r["raw"][lab].get("leaf")]
    return cnt, n, dec_t, all_t, leaf


def ratio_vs_py(rows, lab):
    """Per-case time ratio lab / goursat.py over the cases both decided (verified answer or
    certificate on each side): (median ratio, sum lab / sum goursat.py, count)."""
    pairs = [(r["raw"][lab]["time"], r["raw"][PYLAB]["time"]) for r in rows
             if r["runs"].get(lab) in DECIDED and r["runs"].get(PYLAB) in DECIDED
             and r["raw"][lab].get("time") and r["raw"][PYLAB].get("time")]
    if not pairs:
        return None
    return (statistics.median(m / max(p, 1e-4) for m, p in pairs),
            sum(m for m, _ in pairs) / max(sum(p for _, p in pairs), 1e-4), len(pairs))


def fmt_t(ts, f=max):
    return f"{f(ts):.3g}" if ts else "-"


def fmt_ratio(rt):
    return "-" if rt is None else f"{rt[0]:.3g} ({rt[1]:.3g}, {rt[2]})"


def write_report(path, saved, date):
    rows, labels = saved["rows"], saved["labels"]
    ld = saved.get("load", {})
    corr = saved.get("mode") == "correctness"
    L = [f"# 265-goursat-stress ({date})" + (": correctness pass; timings pending a quiet machine" if corr else ""), "",
         *(["**Correctness pass; timings pending a quiet machine.** This run asks only whether an answer is "
            "right: wrong answers and verdicts do not depend on machine load, times and timeouts do. The "
            "times below are NOT a baseline, and every timeout is provisional (it may close on a quiet "
            "machine, or within the timed run's own limit). The timed baseline is "
            "`python3 tools/goursat/goursat_stress.py` on a quiet machine.", ""] if corr else []),
         f"{len(rows)} cases. {saved.get('mathilda', '?')}; per-case timeout {saved['timeout']} s "
         f"for Mathilda and sympy, {saved['jobs']} parallel processes; wall {saved['wall'] / 60:.1f} min.",
         f"Machine load (1-min average, {ld.get('cpus', '?')} CPUs): start {ld.get('start')}, mean "
         f"{ld.get('mean')}, max {ld.get('max')}, end {ld.get('end')} -- other agents were building "
         "and testing on the same machine, so times are noisy upper bounds."
         + (" The src/internal .m stamp changed during the run (a package was edited mid-run)."
            if saved.get("stamp", {}).get("start") != saved.get("stamp", {}).get("end") else ""), "",
         "Every Mathilda and sympy answer is judged by PossibleZeroQ[D[F, x] - f] in Mathilda: zero is "
         "`verified`, nonzero is `wrong`, anything else is `unverified`. `cert`: certified "
         "non-elementary. goursat.py "
         + ("is read from tests/data/generalised_goursat/sympy_*.tsv (bench_sympy.py, 120 s limit, its own "
            "numerical verification; no answer text, so not PossibleZeroQ-judged)."
            if saved.get("goursatpy") == "tsv" else "was re-run through adapters.run_goursatpy and "
            "judged like the others."), "",
         "Times: `median s` over decided runs (verified answer or certificate), `max s` over all runs "
         "(a timeout counts its full limit). `x goursat.py`: median per-case time ratio method / "
         "goursat.py over the cases both decided (in brackets: total-time ratio, and the number of "
         "such cases); below 1 means faster than goursat.py.", ""]
    if saved.get("notes"):              # hand-written findings stored in results/<run>.json ("notes")
        L += [saved["notes"].rstrip(), ""]
    hdr = ("| run | n | verified | wrong | unverified | cert | false cert | uneval | timeout | error | "
           "median s | max s | median leaves | x goursat.py |")
    sep = "|---" * 14 + "|"

    def line(lab, rs, first):
        cnt, n, dec_t, all_t, leaf = stats_of(rs, lab)
        if not n:
            return None
        fc = sum(1 for r in rs if r["runs"].get(lab) == "CERT" and r["expect"] == "elem")
        rt = ratio_vs_py(rs, lab) if lab != PYLAB else None
        ver = f"{cnt['ELEM'] + cnt['SPECIAL']}" + (f" ({cnt['SPECIAL']} special)" if cnt["SPECIAL"] else "")
        return (f"| {first}{lab} | {n} | {ver} | {cnt['WRONG']} | "
                f"{cnt['UNVERIFIED']} | {cnt['CERT']} | {fc} | {cnt['UNEVAL']} | {cnt['TIMEOUT']} | "
                f"{cnt['ERROR']} | {fmt_t(dec_t, statistics.median)} | {fmt_t(all_t)} | "
                f"{fmt_t(leaf, statistics.median)} | {fmt_ratio(rt)} |")

    L += ["## All cases", "", hdr, sep]
    L += [x for x in (line(lab, rows, "") for lab in labels) if x]
    if saved.get("goursatpy") == "tsv":
        L += ["", "goursat.py `verified` = its own numerical check passed (defect <= 1e-20).", ""]

    L += ["## By family", "", "Families pool the corpora they occur in (the stress families pool "
          "stress_2026, stress_7 and stress_11).", "",
          "| family | expect | " + hdr[2:], "|---" * 16 + "|"]
    fams = collections.OrderedDict()
    for r in rows:
        fams.setdefault(r["family"], []).append(r)
    for fam, rs in fams.items():
        ex = collections.Counter(r["expect"] for r in rs)
        exs = " ".join(f"{k} {v}" for k, v in ex.most_common())
        for lab in labels:
            ln = line(lab, rs, "")
            if ln:
                L.append(f"| {fam} | {exs} " + ln)
    L += ["", "## By corpus (GeneralisedGoursat vs goursat.py)", "",
          "| corpus | n | GG verified | GG cert | GG timeout | goursat.py elementary | goursat.py "
          "non-elementary | goursat.py timeout | GG median s | goursat.py median s | x goursat.py |",
          "|---" * 11 + "|"]
    gg = "Mathilda:GeneralisedGoursat"
    for corp in CORPORA:
        rs = [r for r in rows if r["corpus"] == corp]
        if not rs:
            continue
        c1, _n, d1, _a, _l = stats_of(rs, gg)
        c2, _n2, d2, _a2, _l2 = stats_of(rs, PYLAB)
        L.append(f"| {corp} | {len(rs)} | {c1['ELEM'] + c1['SPECIAL']} | {c1['CERT']} | {c1['TIMEOUT']} | "
                 f"{c2['ELEM']} | {c2['CERT']} | {c2['TIMEOUT']} | {fmt_t(d1, statistics.median)} | "
                 f"{fmt_t(d2, statistics.median)} | {fmt_ratio(ratio_vs_py(rs, gg))} |")

    # hard failures, explicitly
    wrong = [(r, lab) for r in rows for lab, cls in r["runs"].items() if lab.startswith("Mathilda:") and cls == "WRONG"]
    L += ["", f"## Wrong answers from Mathilda ({len(wrong)})", "",
          "An answer PossibleZeroQ[D[F, x] - f] rejects (or a garbage answer carrying $Failed or a "
          "package-private symbol). Each is a wrong answer.", ""]
    for r, lab in wrong:
        raw = r["raw"][lab]
        L.append(f"- `{r['id']}` ({r['label']}) {lab}, {raw.get('time') or 0:.2f} s: "
                 f"`{r['f'][:200]}` -> `{(raw.get('text') or '')[:240]}`")
    fc = [(r, lab) for r in rows for lab, cls in r["runs"].items()
          if lab.startswith("Mathilda:") and cls == "CERT" and r["expect"] == "elem"]
    L += ["", f"## False certificates ({len(fc)})", "",
          "A non-elementary certificate on a case the corpus constructs as elementary.", ""]
    for r, lab in fc:
        L.append(f"- `{r['id']}` ({r['label']}) {lab}: `{r['f'][:200]}`")
    unv = [(r, lab) for r in rows for lab, cls in r["runs"].items() if lab.startswith("Mathilda:") and cls == "UNVERIFIED"]
    L += ["", f"## Unverified Mathilda answers ({len(unv)})", "",
          "PossibleZeroQ did not decide (timed out, or D[F, x] - f does not evaluate numerically).", ""]
    for r, lab in unv[:40]:
        L.append(f"- `{r['id']}` {lab} ({r['raw'][lab].get('pzq')}): `{r['f'][:160]}`")

    # verdicts against the corpus labels (and against goursat.py)
    dis = []
    for r in rows:
        for lab, cls in r["runs"].items():
            if r["expect"] == "elem" and cls == "CERT":
                dis.append((r, lab, "certified non-elementary, corpus says elementary"))
            elif r["expect"] == "nonelem" and cls == "ELEM":
                dis.append((r, lab, "verified ELEMENTARY answer, corpus says non-elementary"))
            elif (r["expect"] == "nonelem" and cls == "WRONG" and lab.startswith("Mathilda:")
                  and not T.special_q("mathilda", r["raw"][lab].get("text") or "")):
                dis.append((r, lab, "ELEMENTARY answer (rejected by PossibleZeroQ), corpus says non-elementary"))
            # a special-function antiderivative of a non-elementary integral is consistent
    L += ["", f"## Verdicts against the corpus labels ({len(dis)})", "",
          "A certificate on an `elem` case, or a verified elementary answer on a `nonelem` case "
          "(a special-function answer to a `nonelem` case is consistent and not listed). For goursat.py "
          "`verified` is its own numerical check.", ""]
    for r, lab, why in dis:
        L.append(f"- `{r['id']}` ({r['label'][:70]}) {lab}: {why}")
    noncert = collections.Counter()
    for r in rows:
        if r["expect"] == "nonelem":
            for lab, cls in r["runs"].items():
                noncert[(lab, cls)] += 1
    if noncert:
        L += ["", "How each run ends on the `nonelem` cases (a certificate or a special-function answer "
              "is a decision; uneval/timeout is not):", "",
              "| run | " + " | ".join(("CERT", "SPECIAL", "ELEM", "UNEVAL", "TIMEOUT", "UNVERIFIED", "WRONG", "ERROR"))
              + " |", "|---" * 9 + "|"]
        for lab in labels:
            L.append(f"| {lab} | " + " | ".join(str(noncert[(lab, k)]) for k in
                     ("CERT", "SPECIAL", "ELEM", "UNEVAL", "TIMEOUT", "UNVERIFIED", "WRONG", "ERROR")) + " |")
    gl = [r for r in rows if r["runs"].get(gg) == "CERT" and r["runs"].get(PYLAB) == "ELEM"]
    gl += [r for r in rows if r["runs"].get(gg) == "ELEM" and r["runs"].get(PYLAB) == "CERT"]
    L += ["", f"GeneralisedGoursat and goursat.py reach opposite verdicts on {len(gl)} cases"
          + (": " + ", ".join(f"`{r['id']}`" for r in gl[:40]) if gl else "") + ".", ""]
    tos = collections.Counter(lab for r in rows for lab, cls in r["runs"].items()
                              if cls == "TIMEOUT" and lab.startswith("Mathilda:"))
    if corr and tos:
        L += [f"Provisional timeouts (at {saved['timeout']} s, under load): "
              + ", ".join(f"{k} {v}" for k, v in tos.items()) + ".", ""]

    # worst offenders
    L += ["", "## Worst offenders" + (" (times under load: indicative only)" if corr else ""), ""]
    for lab in labels:
        if not lab.startswith("Mathilda:"):
            continue
        slow = sorted((r for r in rows if r["runs"].get(lab) in DECIDED and r["raw"][lab].get("time")),
                      key=lambda r: -r["raw"][lab]["time"])[:10]
        L += [f"### {lab}: slowest decided", "", "| case | label | s | class | goursat.py s |", "|---|---|---|---|---|"]
        for r in slow:
            py = r["raw"].get(PYLAB, {})
            L.append(f"| `{r['id']}` | {r['label'][:60]} | {r['raw'][lab]['time']:.2f} | {r['runs'][lab]} | "
                     f"{py.get('time', 0) or 0:.2f} ({r['runs'].get(PYLAB, '-')}) |")
        lost = [r for r in rows if r["runs"].get(PYLAB) in DECIDED and r["runs"].get(lab) not in DECIDED
                and lab in r["runs"]]
        L += ["", f"{lab}: {len(lost)} cases goursat.py decides and this run does not "
              f"({', '.join(f'{k} {v}' for k, v in collections.Counter(r['runs'][lab] for r in lost).most_common())}).", ""]
        if lab == gg:
            for r in lost[:80]:
                py = r["raw"][PYLAB]
                L.append(f"- `{r['id']}` ({r['label'][:60]}): {r['runs'][lab]} "
                         f"{r['raw'][lab].get('time') or 0:.1f} s; goursat.py {r['runs'][PYLAB]} {py.get('time') or 0:.1f} s")
            L.append("")
        slower = []
        for r in rows:
            if r["runs"].get(lab) in DECIDED and r["runs"].get(PYLAB) in DECIDED:
                tm, tp = r["raw"][lab].get("time") or 0, r["raw"][PYLAB].get("time") or 0
                if tp > 0 and tm > tp and tm > 1:
                    slower.append((tm / tp, r, tm, tp))
        slower.sort(key=lambda x: -x[0])
        if slower:
            L += [f"{lab}: slower than goursat.py on {len(slower)} decided cases (over 1 s); the worst:", ""]
            for q, r, tm, tp in slower[:10]:
                L.append(f"- `{r['id']}`: {tm:.2f} s vs {tp:.2f} s ({q:.1f}x)")
            L.append("")
    py_only = [r for r in rows if r["runs"].get(gg) in DECIDED and r["runs"].get(PYLAB) not in DECIDED and PYLAB in r["runs"]]
    L += [f"GeneralisedGoursat decides {len(py_only)} cases goursat.py does not "
          f"({', '.join(f'{k} {v}' for k, v in collections.Counter(r['runs'][PYLAB] for r in py_only).most_common())}).", ""]
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


def history_entry(saved):
    rows = saved["rows"]
    h = {"date": saved["date"], "mode": saved.get("mode", "timed"), "mathilda": saved.get("mathilda"),
         "cases": len(rows),
         "timeout": saved["timeout"], "load_mean": saved.get("load", {}).get("mean"),
         "verified": {}, "cert": {}, "wrong": {}, "timeout_runs": {}, "median_s": {}}
    for lab in saved["labels"]:
        cnt, n, dec_t, _a, _l = stats_of(rows, lab)
        h["verified"][lab] = cnt["ELEM"] + cnt["SPECIAL"]
        h["cert"][lab] = cnt["CERT"]
        h["wrong"][lab] = cnt["WRONG"]
        h["timeout_runs"][lab] = cnt["TIMEOUT"]
        h["median_s"][lab] = round(statistics.median(dec_t), 4) if dec_t else None
    return h


if __name__ == "__main__":
    main()
