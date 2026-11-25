#!/usr/bin/env python3
"""shootout.py -- the integration shootout: every applicable Mathilda Integrate
method and four other systems on one corpus (benchmarks/256-integration-shootout).

    python3 tools/integration_shootout/shootout.py                    # everything
    python3 tools/integration_shootout/shootout.py --quick --systems mathilda,fricas --check
    python3 tools/integration_shootout/shootout.py --classes M,T --limit 50 --jobs 6

RUNS.  Each case of corpus.tsv runs through
  * Mathilda: Method -> Automatic and ParallelMixedTower for every class, plus
    RischTranscendental (classes R, T), Trager (A), Bronstein (M, B);
  * FriCAS, Axiom, Maxima, sympy, the Wolfram Language (wolframscript) and Reduce (PSL),
    when installed.
VERIFICATION.  Every answer of every system is judged by one test, in Mathilda:
PossibleZeroQ[D[F, x] - f] (crosscheck.py; answers of the other systems are
translated into Mathilda first): zero -> verified, nonzero -> WRONG, anything
else -> UNVERIFIED.  Each system also checks its own answer at three points
(adapters.py); that check is recorded ("nums") but decides nothing.

PER-RUN CLASS:  ELEM (verified, elementary) | SPECIAL (verified, uses a
non-elementary special function) | WRONG (PossibleZeroQ[D[F, x] - f] is False) |
UNVERIFIED (PossibleZeroQ did not decide, or the answer would not translate) |
CERT (Mathilda proved it non-elementary: Integrate::nonelem) | UNEVAL |
TIMEOUT | ERROR.

CROSS-RUN VERDICTS (Mathilda runs only):
  FALSE_CERT  CERT on a case whose expect is "elem"              [hard fail]
  DISAGREE    CERT while another run, any method or system, has a verified
              ELEM answer                                         [hard fail]
  WRONG       a Mathilda answer PossibleZeroQ[D[F, x] - f] rejects [hard fail]
  AUTO-GAP    a dedicated method closes, Automatic does not       [reported]
  ONLY-<sys>  exactly one system closes the case (ELEM)          [reported]
Another system's WRONG is listed as that system's bug, never as ours.

--check exits 1 on a hard failure, or when a closed count (ELEM + SPECIAL per
Mathilda method and class) falls below baseline.json (a ratchet;
--update-baseline rewrites it from this run).

CACHE.  Every finished run is appended to results/cache.jsonl keyed by system,
method, case id and integrand (Mathilda runs also by the newest mtime of the
binary and the src/internal .m packages), so an
interrupted run resumes and an unchanged external system is not re-run;
--fresh ignores the cache.
"""
import argparse
import collections
import concurrent.futures as cf
import datetime
import hashlib
import json
import os
import re
import statistics
import sys
import tempfile
import threading

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import adapters as A  # noqa: E402
import translate as T
from easy_classify import easy_ids  # noqa: E402
import crosscheck as X  # noqa: E402

BENCH = os.path.join(A.REPO, "benchmarks", "256-integration-shootout")
CORPUS = os.path.join(BENCH, "corpus.tsv")
EXTERNAL = ("fricas", "axiom", "maxima", "sympy", "wolfram", "reduce")
LABEL = {"fricas": "FriCAS", "axiom": "Axiom", "maxima": "Maxima", "sympy": "sympy", "wolfram": "Wolfram", "reduce": "Reduce"}
METHODS = {"R": ["Automatic", "ParallelMixedTower", "RischTranscendental"],
           "T": ["Automatic", "ParallelMixedTower", "RischTranscendental"],
           "A": ["Automatic", "ParallelMixedTower", "Trager"],
           "M": ["Automatic", "ParallelMixedTower", "Bronstein"],
           "B": ["Automatic", "ParallelMixedTower", "Bronstein"]}
ALL_METHODS = ["Automatic", "ParallelMixedTower", "RischTranscendental", "Trager", "Bronstein"]
CLOSED = ("ELEM", "SPECIAL")



def mathilda_stamp():
    """The version of Mathilda a cached run belongs to: the binary's mtime AND
    the newest mtime under src/internal (the .m packages -- Bronstein, Trager,
    ParallelMixedSpecial -- are loaded at run time, so editing one changes the
    answers without touching the binary)."""
    newest = os.path.getmtime(A.MATHILDA)
    for root, _dirs, files in os.walk(os.path.join(A.REPO, "src", "internal")):
        for fn in files:
            if fn.endswith(".m"):
                newest = max(newest, os.path.getmtime(os.path.join(root, fn)))
    return str(int(newest))

def load_corpus(path=None):
    path = path or CORPUS
    cases = []
    for ln in open(path):
        p = ln.rstrip("\n").split("\t")
        if len(p) < 10 or p[0] != "CASE":
            continue
        pts = p[7].strip("{}")
        c = {"id": p[1], "class": p[2], "family": p[3], "expect": p[4], "f": p[5],
             "F": None if p[6] == "None" else p[6], "points": p[7],
             "points_list": [s.strip() for s in pts.split(",")], "fullform": p[9]}
        c["tree"] = T.parse(c["fullform"])
        c["key"] = hashlib.sha1(c["f"].encode()).hexdigest()[:12]
        cases.append(c)
    return cases


def run_label(system, method):
    return f"Mathilda:{method}" if system == "mathilda" else LABEL[system]


def classify_run(system, raw, xv=None):
    """xv: the PossibleZeroQ[D[F, x] - f] verdict for the answer, which decides it; None
    only under --no-crosscheck, when the system's own point check is used instead."""
    st = raw["status"]
    if st == "UNEVAL":
        return "CERT" if raw.get("cert") else "UNEVAL"
    if st != "ANSWER":
        return st
    nums = raw.get("nums") or []
    if xv is not None:
        if xv == "nonzero":
            return "WRONG"
        if xv != "zero":
            return "UNVERIFIED"
        nums = ["ok"]
    if "bad" in nums:
        return "WRONG"
    if "ok" not in nums:
        return "UNVERIFIED"
    return "SPECIAL" if T.special_q("mathilda" if system == "mathilda" else system, raw.get("text", "")) else "ELEM"


def verdicts(case, runs):
    """runs: {label: cls}.  Returns the list of cross-run verdict strings."""
    out = []
    elem_by = [lab for lab, cls in runs.items() if cls == "ELEM"]
    for lab, cls in runs.items():
        if not lab.startswith("Mathilda:"):
            continue
        if cls == "WRONG":
            out.append(f"WRONG {lab}")
        if cls == "CERT":
            if case["expect"] == "elem":
                out.append(f"FALSE_CERT {lab}")
            elif elem_by:
                out.append(f"DISAGREE {lab} vs {elem_by[0]}")
        if cls in CLOSED and case["expect"] == "nonelem" and cls == "ELEM":
            out.append(f"CORPUS? {lab} closed a theorem-nonelem case")
    auto = runs.get("Mathilda:Automatic")
    if auto is not None and auto not in CLOSED:
        for lab, cls in runs.items():
            if lab.startswith("Mathilda:") and lab != "Mathilda:Automatic" and cls in CLOSED:
                out.append(f"AUTO-GAP {lab}")
                break
    sys_elem = set()
    for lab, cls in runs.items():
        if cls == "ELEM":
            sys_elem.add("Mathilda" if lab.startswith("Mathilda:") else lab)
    if len(sys_elem) == 1 and len({("Mathilda" if l.startswith("Mathilda:") else l) for l in runs}) > 1:
        out.append(f"ONLY-{next(iter(sys_elem))}")
    return out


HARD = ("WRONG", "FALSE_CERT", "DISAGREE")


class Cache:
    def __init__(self, path, fresh):
        self.path, self.lock, self.d = path, threading.Lock(), {}
        if not fresh and os.path.exists(path):
            for ln in open(path):
                try:
                    r = json.loads(ln)
                    self.d[r["ckey"]] = r
                except (ValueError, KeyError):
                    pass

    def get(self, k):
        return self.d.get(k)

    def put(self, k, rec):
        rec = dict(rec, ckey=k)
        with self.lock:
            self.d[k] = rec
            with open(self.path, "a") as fh:
                fh.write(json.dumps(rec) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", default="mathilda," + ",".join(EXTERNAL))
    ap.add_argument("--methods", default=",".join(ALL_METHODS))
    ap.add_argument("--classes", default="R,T,A,M,B")
    ap.add_argument("--families", default=None, help="comma-separated family prefixes")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--quick", action="store_true", help="every 8th case (the nightly subset)")
    ap.add_argument("--timeout", type=int, default=30)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--update-baseline", action="store_true")
    ap.add_argument("--fresh", action="store_true")
    ap.add_argument("--include-easy", action="store_true",
                    help="run the cases <bench>/easy.tsv classes as easy too (rational, rationalisable by a "
                    "substitution, closed by DerivativeDivides); they are left out by default")
    ap.add_argument("--no-crosscheck", action="store_true",
                    help="judge answers by their own system's check only (no PossibleZeroQ cross-check)")
    ap.add_argument("--bench", default=None, help="the benchmark folder (corpus.tsv, results/, REPORT.md, "
                    "baseline.json); default benchmarks/256-integration-shootout")
    ap.add_argument("--report", default=None, help="REPORT.md path (default: the benchmark folder; "
                    "--quick/--limit runs write REPORT-partial.md)")
    ap.add_argument("--retabulate", action="store_true",
                    help="rewrite REPORT.md from the latest results/<date>.json without running anything, "
                    "leaving out the easy cases (unless --include-easy)")
    a = ap.parse_args()
    global BENCH, CORPUS
    if a.bench:
        BENCH = os.path.abspath(a.bench)
        CORPUS = os.path.join(BENCH, "corpus.tsv")
    if a.retabulate:
        retabulate(a)
        return

    if not A.available("mathilda"):
        print("SKIP: ./Mathilda not built")
        sys.exit(77)
    systems = [s for s in a.systems.split(",") if s]
    missing = [s for s in systems if not A.available(s)]
    for s in missing:
        print(f"note: {s} not available, skipped")
    systems = [s for s in systems if s not in missing]
    methods = set(a.methods.split(","))
    classes = set(a.classes.split(","))

    cases = [c for c in load_corpus() if c["class"] in classes]
    if not a.include_easy:
        easy = easy_ids(BENCH)
        if easy:
            print(f"note: {sum(c['id'] in easy for c in cases)} easy cases left out (easy.tsv; --include-easy keeps them)")
        cases = [c for c in cases if c["id"] not in easy]
    if a.families:
        fams = a.families.split(",")
        cases = [c for c in cases if any(c["family"].startswith(f) for f in fams)]
    if a.quick:
        cases = cases[::8]
    if a.limit:
        cases = cases[:a.limit]
    os.makedirs(os.path.join(BENCH, "results"), exist_ok=True)
    cache = Cache(os.path.join(BENCH, "results", "cache.jsonl"), a.fresh)
    mstamp = mathilda_stamp()

    # work units: (system, method, [cases]) chunks
    units = []
    for s in systems:
        if s == "mathilda":
            for m in ALL_METHODS:
                if m not in methods:
                    continue
                cs = [c for c in cases if m in METHODS[c["class"]]]
                units += [("mathilda", m, cs[i:i + 8]) for i in range(0, len(cs), 8)]
        else:
            step = 20 if s in ("wolfram", "sympy") else 4
            units += [(s, None, cases[i:i + step]) for i in range(0, len(cases), step)]

    def ckey(s, m, c):
        return "|".join([s, m or "", c["id"], c["key"], mstamp if s == "mathilda" else ""])

    wolf_sem = threading.Semaphore(2)
    tmpdir = tempfile.mkdtemp(prefix="shootout-")
    done = collections.Counter()
    total = sum(len(u[2]) for u in units)
    plock = threading.Lock()

    def work(u):
        s, m, cs = u
        todo = [c for c in cs if cache.get(ckey(s, m, c)) is None]
        if todo:
            if s == "wolfram":
                with wolf_sem:
                    got = A.run_system(s, todo, a.timeout, tmpdir, m)
            else:
                got = A.run_system(s, todo, a.timeout, tmpdir, m)
            for c in todo:
                r = got.get(c["id"]) or {"id": c["id"], "status": "ERROR", "time": None, "nums": [],
                                          "cert": False, "text": "", "detail": "no record"}
                cache.put(ckey(s, m, c), r)
        with plock:
            done["n"] += len(cs)
            if done["n"] % 100 < len(cs):
                print(f"  {done['n']}/{total} runs", flush=True)

    t0 = datetime.datetime.now()
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        list(ex.map(work, units))
    xc = {} if a.no_crosscheck else crosscheck(cases, systems, methods, cache, ckey, mstamp, tmpdir, a.jobs)
    wall = (datetime.datetime.now() - t0).total_seconds()

    # assemble
    rows = []
    for c in cases:
        runs, raw = {}, {}
        for s in systems:
            ms = [m for m in ALL_METHODS if m in methods and m in METHODS[c["class"]]] if s == "mathilda" else [None]
            for m in ms:
                r = cache.get(ckey(s, m, c))
                if r is None:
                    continue
                lab = run_label(s, m)
                xv = xc.get(ckey(s, m, c))
                runs[lab] = classify_run(s, r, xv)
                raw[lab] = {k: r.get(k) for k in ("status", "time", "nums", "cert", "detail")}
                if xv is not None:
                    raw[lab]["pzq"] = xv
                raw[lab]["text"] = (r.get("text") or "")[:300]
        rows.append({"id": c["id"], "class": c["class"], "family": c["family"], "expect": c["expect"],
                     "f": c["f"], "runs": runs, "raw": raw, "verdicts": verdicts(c, runs)})

    full = not (a.quick or a.limit or a.families or a.classes != "R,T,A,M,B")
    report = a.report or os.path.join(BENCH, "REPORT.md" if full else "REPORT-partial.md")
    labels = [run_label("mathilda", m) for m in ALL_METHODS if m in methods and "mathilda" in systems] + \
             [LABEL[s] for s in systems if s != "mathilda"]
    write_report(report, rows, labels, wall, a)
    stamp = datetime.date.today().isoformat()
    if full:
        versions = {(run_label("mathilda", "Automatic").split(":")[0] if sy == "mathilda" else LABEL[sy]):
                    A.version(sy) for sy in systems}
        json.dump({"date": stamp, "wall": wall, "systems": versions, "rows": rows},
                  open(os.path.join(BENCH, "results", f"{stamp}.json"), "w"), indent=0)
        with open(os.path.join(BENCH, "history.jsonl"), "a") as fh:
            fh.write(json.dumps({"date": stamp, "cases": len(rows), "closed": closed_counts(rows, labels)}) + "\n")

    hard = [(r["id"], v) for r in rows for v in r["verdicts"] if v.split()[0] in HARD]
    for cid, v in hard:
        print(f"*** {v}  [{cid}]")
    fail = bool(hard)
    bpath = os.path.join(BENCH, "baseline.json")
    counts = mathilda_closed(rows)
    if a.update_baseline:
        base = json.load(open(bpath)) if os.path.exists(bpath) else {}
        base[mode_key(a)] = counts
        json.dump(base, open(bpath, "w"), indent=1, sort_keys=True)
    elif a.check and os.path.exists(bpath):
        base = json.load(open(bpath)).get(mode_key(a), {})
        for k, v in base.items():
            if counts.get(k, 0) < v:
                print(f"*** RATCHET {k}: {counts.get(k, 0)} closed < baseline {v}")
                fail = True
    print(f"report: {report}  ({len(rows)} cases, {wall:.0f} s)")
    sys.exit(1 if (a.check and fail) else 0)


def crosscheck(cases, systems, methods, cache, ckey, mstamp, tmpdir, jobs):
    """Judge, in Mathilda by PossibleZeroQ, every answer of every system.
    Returns {run ckey: verdict}.  Verdicts are cached in the run cache under
    crosscheck.ckey (keyed by the Mathilda stamp and the answer text)."""
    todo, out = [], {}
    for c in cases:
        for s in systems:
            ms = [m for m in ALL_METHODS if m in methods and m in METHODS[c["class"]]] if s == "mathilda" else [None]
            for m in ms:
                k = ckey(s, m, c)
                r = cache.get(k)
                if r is None or r.get("status") != "ANSWER":
                    continue
                if s == "mathilda" and r.get("nums") == ["bad"]:
                    continue        # garbage ($Failed / a Private` symbol): wrong as it stands
                xk = X.ckey(k + "|" + mstamp, r.get("text") or "")
                hit = cache.get(xk)
                if hit is not None:
                    out[k] = hit["verdict"]
                    continue
                try:
                    wl = r["text"] if s in ("mathilda", "wolfram") else X.to_wl(s, r["text"])
                except X.Untranslatable as e:
                    cache.put(xk, {"id": c["id"], "verdict": "untranslatable", "detail": str(e)})
                    out[k] = "untranslatable"
                    continue
                todo.append((k, xk, c, wl))
    if not todo:
        return out
    print(f"  cross-check: {len(todo)} answers to judge in Mathilda", flush=True)
    chunks = [todo[i:i + 10] for i in range(0, len(todo), 10)]

    def work(chunk):
        # ids inside a batch must be unique: number the items
        items = [(f"x{j}", t[2]["f"], t[3]) for j, t in enumerate(chunk)]
        got = X.run_batch(items, tmpdir)
        for j, (k, xk, c, wl) in enumerate(chunk):
            v = got.get(f"x{j}", "undecided")
            cache.put(xk, {"id": c["id"], "verdict": v, "detail": ""})
            out[k] = v

    with cf.ThreadPoolExecutor(max_workers=jobs) as ex:
        list(ex.map(work, chunks))
    return out


def mode_key(a):
    return "quick" if a.quick else "full"


def mathilda_closed(rows):
    c = collections.Counter()
    for r in rows:
        for lab, cls in r["runs"].items():
            if lab.startswith("Mathilda:") and cls in CLOSED:
                c[f"{lab}/{r['class']}"] += 1
    return dict(c)


def closed_counts(rows, labels):
    return {lab: sum(1 for r in rows if r["runs"].get(lab) in CLOSED) for lab in labels}


def pct(n, d):
    return f"{100.0 * n / d:.0f}%" if d else "-"


def retabulate(a):
    """REPORT.md from the latest saved run: the same tables, without the easy cases"""
    rdir = os.path.join(BENCH, "results")
    date = sorted(f[:-5] for f in os.listdir(rdir) if re.match(r"\d{4}-\d\d-\d\d\.json$", f))[-1]
    saved = json.load(open(os.path.join(rdir, date + ".json")))
    easy = set() if a.include_easy else easy_ids(BENCH)
    rows = [r for r in saved["rows"] if r["id"] not in easy]
    labels = []
    for r in saved["rows"]:
        labels += [k for k in r["runs"] if k not in labels]
    mathilda = [run_label("mathilda", m) for m in ALL_METHODS]
    labels = [l for l in mathilda if l in labels] + [l for l in labels if l not in mathilda]
    path = a.report or os.path.join(BENCH, "REPORT.md")
    write_report(path, rows, labels, saved["wall"], a, run_date=date,
                 note=f"Re-tabulated from the run of {date}: {len(saved['rows']) - len(rows)} easy cases "
                      "(easy.tsv: rational, rationalised by a substitution, or closed by DerivativeDivides) "
                      "are left out." if easy else "")
    print(f"{path}: {len(rows)} of {len(saved['rows'])} cases (run {date})")


def write_report(path, rows, labels, wall, a, run_date=None, note=""):
    L = [f"# Integration shootout: {os.path.basename(BENCH)} ({run_date or datetime.date.today().isoformat()})", "",
         *([note, ""] if note else []),
         f"{len(rows)} cases, timeout {a.timeout} s per run, wall {wall / 60:.1f} min. "
         "Every answer of every system is judged by PossibleZeroQ[D[F, x] - f] in Mathilda (the other "
         "systems' answers translated into Mathilda first): zero is verified, nonzero is `wrong`, "
         "anything else -- undecided, untranslatable, or an answer text cut when the run was recorded -- "
         "is `unverified`. `elem`/`special`: verified elementary / special-function answer; `cert`: "
         "Mathilda proved non-elementary.", ""]
    classes = [c for c in "RTAMB" if any(r["class"] == c for r in rows)]
    name = {"R": "rational", "T": "transcendental", "A": "algebraic", "M": "mixed (radical on top)",
            "B": "radical at the bottom"}
    L += ["## Closed (elementary + special) by class", "",
          "| run | " + " | ".join(f"{c} {name[c]} ({sum(1 for r in rows if r['class'] == c)})" for c in classes)
          + " | all |", "|---" * (len(classes) + 2) + "|"]
    for lab in labels:
        cells = []
        tot_n = tot_d = 0
        for c in classes:
            rs = [r for r in rows if r["class"] == c and lab in r["runs"]]
            n = sum(1 for r in rs if r["runs"][lab] in CLOSED)
            cells.append(f"{pct(n, len(rs))} ({n})" if rs else "·")
            tot_n += n
            tot_d += len(rs)
        L.append(f"| {lab} | " + " | ".join(cells) + f" | {pct(tot_n, tot_d)} ({tot_n}) |")
    L += ["", "## Outcome counts", "",
          "| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for lab in labels:
        cnt = collections.Counter(r["runs"][lab] for r in rows if lab in r["runs"])
        ts = [r["raw"][lab]["time"] for r in rows if r["runs"].get(lab) in CLOSED and r["raw"][lab].get("time")]
        med = f"{statistics.median(ts):.2f}" if ts else "-"
        L.append(f"| {lab} | {cnt['ELEM']} | {cnt['SPECIAL']} | {cnt['CERT']} | {cnt['UNEVAL']} | "
                 f"{cnt['TIMEOUT']} | {cnt['ERROR']} | {cnt['UNVERIFIED']} | {cnt['WRONG']} | {med} |")
    xrows = [(lab, r["raw"][lab]) for r in rows for lab in labels
             if lab in r["raw"] and r["raw"][lab].get("pzq") is not None]
    if xrows:
        L += ["", "## PossibleZeroQ verdicts", "",
              "Every answer, by PossibleZeroQ[D[F, x] - f] in Mathilda. `own ok`: of the verified answers, "
              "those the system's own three-point check also passed; `own disagrees`: answers where the "
              "own check and PossibleZeroQ disagree (one passes, the other does not).", "",
              "| run | answers | zero (verified) | nonzero (wrong) | undecided | untranslatable | cut | "
              "own ok | own disagrees |",
              "|---|---|---|---|---|---|---|---|---|"]
        for lab in labels:
            xs = [x for l2, x in xrows if l2 == lab]
            if not xs:
                continue
            cc = collections.Counter(x["pzq"] if x["pzq"] in ("zero", "nonzero", "untranslatable", "cut")
                                     else "undecided" for x in xs)

            def own_ok(x):
                nums = x.get("nums") or []
                return bool(nums) and all(v == "ok" for v in nums)
            ok = sum(1 for x in xs if x["pzq"] == "zero" and own_ok(x))
            dis = sum(1 for x in xs if (x["pzq"] == "zero") != own_ok(x) and x["pzq"] in ("zero", "nonzero"))
            L.append(f"| {lab} | {len(xs)} | {cc['zero']} | {cc['nonzero']} | {cc['undecided']} | "
                     f"{cc['untranslatable']} | {cc['cut']} | {ok} | {dis} |")
    # by family, closed for the best Mathilda run vs each external
    L += ["", "## Closed by family (Mathilda = best of its methods)", "",
          "| class | family | n | Mathilda | " + " | ".join(l for l in labels if not l.startswith("Mathilda:")) + " |",
          "|---" * (4 + sum(1 for l in labels if not l.startswith("Mathilda:"))) + "|"]
    fams = collections.OrderedDict()
    for r in rows:
        fams.setdefault((r["class"], r["family"]), []).append(r)
    for (c, f), rs in fams.items():
        best = sum(1 for r in rs if any(cls in CLOSED for lab, cls in r["runs"].items() if lab.startswith("Mathilda:")))
        ext = [str(sum(1 for r in rs if r["runs"].get(l) in CLOSED)) for l in labels if not l.startswith("Mathilda:")]
        L.append(f"| {c} | {f} | {len(rs)} | {best} | " + " | ".join(ext) + " |")
    ver = collections.defaultdict(list)
    for r in rows:
        for v in r["verdicts"]:
            ver[v.split()[0]].append((r, v))
    for kind in ("WRONG", "FALSE_CERT", "DISAGREE", "CORPUS?", "AUTO-GAP"):
        L += ["", f"## {kind} ({len(ver.get(kind, []))})", ""]
        for r, v in ver.get(kind, [])[:60]:
            L.append(f"- `{r['id']}` {v}: `{r['f'][:220]}`")
    L += ["", "## Only one system closes it", ""]
    only = collections.Counter(v for r in rows for v in r["verdicts"] if v.startswith("ONLY-"))
    for k, n in only.most_common():
        L.append(f"- {k[5:]}: {n}")
    L += ["", "## Wrong answers from other systems (their bugs)", ""]
    for r in rows:
        for lab, cls in r["runs"].items():
            if cls == "WRONG" and not lab.startswith("Mathilda:"):
                L.append(f"- {lab} `{r['id']}`: `{r['f'][:160]}`")
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
