#!/usr/bin/env python3
"""compare_fricas.py -- Mathilda's decision procedures against FriCAS's integrate,
read off the shootout's saved runs (results/<date>.json of benchmarks 256, 257, 259).

    python3 tools/integration_shootout/compare_fricas.py [--date 2026-11-20] [--no-progress] \
        [--out benchmarks/BRONSTEIN_VS_FRICAS.md]

The shootout report counts CLOSED answers only.  For a decision procedure the
question is also how often it DECIDES, so this report counts both sides' "no":
  Mathilda  closed = ELEM + SPECIAL (verified); no = CERT (Integrate::nonelem, a
            certificate); undecided = UNEVAL, TIMEOUT, ERROR.
  FriCAS    closed = a verified answer; no = integrate(...) returned with an
            unevaluated integral(...) -- FriCAS's Risch implementation does that
            only when it has proved the integral non-elementary; undecided = an
            error (mostly "implementation incomplete"), a timeout, or an answer
            it could not verify.  (The shootout's FriCAS UNEVAL pattern also
            matches %A/%B, the dummies of an answer summed over roots; those are
            re-split here from the raw text.)
"Mathilda" is the Bronstein-family method of each class: RischTranscendental
(classes R, T -- Bronstein, Symbolic Integration I), Trager (A), Bronstein (M, B:
the mixed algebraic-transcendental package); "Automatic" is the full cascade.

Soundness: a "no" on a case whose expect is "elem" (f = D[F]) is false; so is a
"no" contradicted by the other system's verified answer.

Speed: on the cases both close, the per-case time each system reports for its
integration (Mathilda's AbsoluteTiming, FriCAS's evaluation time; neither
includes start-up), as medians, the 90th percentile, and the geometric mean of
the per-case ratio Mathilda/FriCAS (< 1: Mathilda faster).
"""
import argparse
import collections
import json
import math
import os
import re
import statistics

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
BENCHES = [("256-integration-shootout", "256: the shootout corpus"),
           ("257-compact-integrals", "257: compact integrands, large integrals"),
           ("259-transcendental-decision", "259: the transcendental decision corpus"),
           ("260-radical-bottom", "260: a radical at the bottom"),
           ("261-compact-integrals-x4", "261: compact integrands, integrals four times larger")]
DEDICATED = {"R": "RischTranscendental", "T": "RischTranscendental", "A": "Trager", "M": "Bronstein", "B": "Bronstein"}
CLASSNAME = {"R": "rational", "T": "transcendental", "A": "algebraic", "M": "mixed, radical on top",
             "B": "radical at the bottom"}
CLOSED = ("ELEM", "SPECIAL")


def mathilda_outcome(row, method):
    cls = row["runs"].get(f"Mathilda:{method}")
    if cls is None:
        return None, None
    t = (row["raw"].get(f"Mathilda:{method}") or {}).get("time")
    if cls in CLOSED:
        return "closed", t
    if cls == "CERT":
        return "no", t
    if cls == "WRONG":
        return "wrong", t
    return "undecided", t


def fricas_outcome(row):
    cls = row["runs"].get("FriCAS")
    raw = row["raw"].get("FriCAS") or {}
    t = raw.get("time")
    if cls is None:
        return None, None
    if cls in CLOSED:
        return "closed", t
    if cls == "WRONG":
        return "wrong", t
    text = raw.get("text") or ""
    if cls == "UNEVAL" and re.search(r"\bintegral\(", text):
        return "no", t
    return "undecided", t


def load(bench, date):
    p = os.path.join(REPO, "benchmarks", bench, "results", f"{date}.json")
    return json.load(open(p))["rows"] if os.path.exists(p) else None


def pct(n, d):
    return f"{100.0 * n / d:.0f}%" if d else "-"


def fmt_t(v):
    if v is None:
        return "-"
    return f"{v * 1000:.0f} ms" if v < 1 else f"{v:.1f} s"


def quant(xs, q):
    xs = sorted(xs)
    if not xs:
        return None
    k = min(len(xs) - 1, max(0, int(math.ceil(q * len(xs))) - 1))
    return xs[k]


# The runs compared, in column order: (label, outcome function).  "Mathilda" is
# the Bronstein-family method of the case's class; ParallelMixedTower is the
# Risch-Norman (parallel) heuristic over a mixed tower, run on every class.
RUNS = [("Mathilda", lambda r: mathilda_outcome(r, DEDICATED[r["class"]])),
        ("PMT", lambda r: mathilda_outcome(r, "ParallelMixedTower")),
        ("Automatic", lambda r: mathilda_outcome(r, "Automatic")),
        ("FriCAS", fricas_outcome)]
MATHILDA_RUNS = ("Mathilda", "PMT", "Automatic")
RUN_LABEL = {"Mathilda": "Mathilda:{dedicated}", "PMT": "Mathilda:ParallelMixedTower",
             "Automatic": "Mathilda:Automatic", "FriCAS": "FriCAS"}


def outcomes(r):
    """{run: outcome} for the runs present on r; a "no" on an elementary-by-
    construction case is "falseno" (it decides nothing)."""
    out = {}
    for name, fn in RUNS:
        o, _ = fn(r)
        if o == "no" and r["expect"] == "elem":
            o = "falseno"
        if o:
            out[name] = o
    return out


def run_class(r, name):
    lab = RUN_LABEL[name].format(dedicated=DEDICATED[r["class"]])
    return r["runs"].get(lab)


def tally(rows, key):
    """key(row) -> group; returns {group: Counter of outcomes per run}"""
    out = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
    for r in rows:
        g = key(r)
        for name, o in outcomes(r).items():
            out[g][name][o] += 1
            out[g][name]["n"] += 1
            # the same SPECIAL / ELEM split for every system (translate.special_q
            # on the system's own answer), so neither side's special-function
            # answers are folded into its elementary count
            if o == "closed":
                out[g][name]["special" if run_class(r, name) == "SPECIAL" else "elem"] += 1
    return out


def table(groups, order, label):
    names = [n for n, _ in RUNS]
    L = [f"| {label} | n | " + " | ".join(f"{s} elementary | {s} special | {s} decided" for s in names) + " |",
         "|---|---|" + "---|---|---|" * len(names)]
    tot = collections.defaultdict(collections.Counter)
    for g in order:
        c = groups[g]
        n = max(c[s]["n"] for s in names)
        cells = []
        for s in names:
            if not c[s]["n"]:
                cells += ["·", "·", "·"]
                continue
            cl = c[s]["closed"]
            el, sp = c[s]["elem"], c[s]["special"]
            de = cl + c[s]["no"]
            cells += [f"{pct(el, n)} ({el})", f"{pct(sp, n)} ({sp})", f"{pct(de, n)} ({de})"]
            tot[s]["closed"] += cl
            tot[s]["elem"] += el
            tot[s]["special"] += sp
            tot[s]["dec"] += de
            tot[s]["present"] = 1
        tot["n"]["n"] += n
        L.append(f"| {g} | {n} | " + " | ".join(cells) + " |")
    n = tot["n"]["n"]
    L.append(f"| **all** | {n} | " + " | ".join(
        (f"**{pct(tot[s][k], n)} ({tot[s][k]})**" if tot[s]["present"] else "·")
        for s in names for k in ("elem", "special", "dec")) + " |")
    return L


def soundness(rows):
    """A "no" means: no ELEMENTARY integral.  It is contradicted only by a
    verified elementary answer (a special-function answer is consistent)."""
    out = []
    for r in rows:
        oc = {}
        for name, fn in RUNS:
            oc[name], _ = fn(r)
        elem = {name: run_class(r, name) == "ELEM" for name, _ in RUNS}
        m_elem = any(elem[n] for n in MATHILDA_RUNS)
        for name in MATHILDA_RUNS + ("FriCAS",):
            who = "Mathilda" if name == "Mathilda" else name
            if oc.get(name) == "no":
                if r["expect"] == "elem":
                    out.append((f"{who} false no (expect elem)", r))
                elif name == "FriCAS" and m_elem:
                    out.append(("FriCAS no vs Mathilda verified elementary answer", r))
                elif name != "FriCAS" and elem["FriCAS"]:
                    out.append((f"{who} no vs FriCAS verified elementary answer", r))
            if oc.get(name) == "wrong":
                out.append((f"{who} answer fails its derivative check", r))
    return out


def speed(rows, name="Mathilda"):
    fn = dict(RUNS)[name]
    pairs = []
    for r in rows:
        m, tm = fn(r)
        f, tf = fricas_outcome(r)
        if m == "closed" and f == "closed" and tm and tf:
            pairs.append((max(tm, 1e-4), max(tf, 1e-3)))
    if not pairs:
        return None
    ms, fs = [p[0] for p in pairs], [p[1] for p in pairs]
    gm = math.exp(statistics.mean(math.log(a / b) for a, b in pairs))
    faster = sum(1 for a, b in pairs if a < b)
    return {"n": len(pairs), "m_med": statistics.median(ms), "f_med": statistics.median(fs),
            "m_p90": quant(ms, 0.9), "f_p90": quant(fs, 0.9), "gm": gm, "m_faster": faster,
            "m_tot": sum(ms), "f_tot": sum(fs), "m_slow": sum(1 for x in ms if x > 1)}


def decided_by(r, name):
    o = outcomes(r).get(name)
    return o in ("closed", "no")


def only_one(rows):
    """cases decided by Mathilda (any method) and not FriCAS, and the reverse,
    by family"""
    c = collections.Counter()
    for r in rows:
        md = any(decided_by(r, n) for n in MATHILDA_RUNS)
        fd = decided_by(r, "FriCAS")
        if md and not fd:
            c[("Mathilda only", r["family"])] += 1
        if fd and not md:
            c[("FriCAS only", r["family"])] += 1
    return c


def method_overlap(rows):
    """within Mathilda: decided by PMT and not the Bronstein-family method, and
    the reverse, by class and family"""
    c = collections.Counter()
    for r in rows:
        pd, dd = decided_by(r, "PMT"), decided_by(r, "Mathilda")
        if pd and not dd:
            c[("PMT only", r["class"], r["family"])] += 1
        if dd and not pd:
            c[("Bronstein-family only", r["class"], r["family"])] += 1
    return c


def fricas_errors(rows):
    c = collections.Counter()
    for r in rows:
        f, _ = fricas_outcome(r)
        if f == "undecided":
            raw = r["raw"].get("FriCAS") or {}
            cls = r["runs"].get("FriCAS")
            det = (raw.get("detail") or "").strip()
            if cls == "TIMEOUT":
                det = "timeout (30 s)"
            elif cls in ("UNEVAL",):
                det = "answer with %A/%B root sums, not verified"
            elif cls == "UNVERIFIED":
                det = "answer not verified"
            if det.startswith("Time:") or "->" in det:
                det = "error, no message (a Lisp-level failure)"
            det = re.sub(r"\s+", " ", det)[:70] or cls
            c[det] += 1
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", default=None, help="results/<date>.json to read (default: newest)")
    ap.add_argument("--no-progress", action="store_true",
                    help="omit the progress table (Mathilda against the corpus's previous saved run)")
    ap.add_argument("--out", default=os.path.join(REPO, "benchmarks", "BRONSTEIN_VS_FRICAS.md"))
    a = ap.parse_args()
    L = []
    allrows = []
    progress = []
    for bench, title in BENCHES:
        rdir = os.path.join(REPO, "benchmarks", bench, "results")
        dates = sorted(f[:-5] for f in os.listdir(rdir) if re.match(r"\d{4}-\d\d-\d\d\.json$", f)) if os.path.isdir(rdir) else []
        date = a.date or (dates[-1] if dates else None)
        rows = load(bench, date) if date else None
        if not rows:
            continue
        rows = [r for r in rows if "FriCAS" in r["runs"]]
        allrows += rows
        L.append(f"## {title} ({len(rows)} cases, run {date})")
        L.append("")
        by_class = tally(rows, lambda r: f"{r['class']} {CLASSNAME[r['class']]}")
        order = [g for g in (f"{k} {CLASSNAME[k]}" for k in "RTAMB") if g in by_class]
        L += table(by_class, order, "class")
        L.append("")
        if bench.startswith("259") or bench.startswith("260"):
            by_fam = tally(rows, lambda r: r["family"])
            L.append("By family (`-elem`: f = D[F], elementary by construction; `-rand`: random):")
            L.append("")
            L += table(by_fam, sorted(by_fam), "family")
            L.append("")
        for nm in ("Mathilda", "PMT"):
            sp = speed(rows, nm)
            if sp:
                L.append(f"Speed, {nm} against FriCAS, on the {sp['n']} cases both close: median "
                         f"{fmt_t(sp['m_med'])} vs {fmt_t(sp['f_med'])}; 90th percentile {fmt_t(sp['m_p90'])} vs "
                         f"{fmt_t(sp['f_p90'])}; total {fmt_t(sp['m_tot'])} vs {fmt_t(sp['f_tot'])}; geometric mean "
                         f"of the ratio {nm}/FriCAS {sp['gm']:.2f}; {nm} faster on {sp['m_faster']} of {sp['n']}; "
                         f"{sp['m_slow']} over 1 s.")
                L.append("")
        prev = [d for d in dates if d < date]
        if not a.no_progress and prev:
            old = load(bench, prev[-1])
            if old:
                def cnt(rs, meth):
                    cl = de = 0
                    for r in rs:
                        o, _ = mathilda_outcome(r, DEDICATED[r["class"]] if meth == "dedicated" else meth)
                        cl += o == "closed"
                        de += o in ("closed", "no")
                    return cl, de
                progress.append((f"{title} ({prev[-1]})", len(old), cnt(old, "dedicated"), cnt(rows, "dedicated"),
                                 cnt(old, "Automatic"), cnt(rows, "Automatic")))
    head = ["# Bronstein's algorithm in Mathilda against FriCAS's integrate", "",
            "Generated by `tools/integration_shootout/compare_fricas.py` from the integration shootout's saved runs; "
            "see its docstring for the definitions. Every answer is verified by the system that gave it "
            "(derivative at three points to 10^-20, or, for FriCAS answers with Ei / li / erf, which it cannot "
            "evaluate numerically, `normalize(D(r, x) - f) = 0`). From 261 on, an answer that fails its own check, "
            "or that its own check cannot evaluate, is translated into Mathilda and judged by `PossibleZeroQ[D[F, x] "
            "- f]`, the test Mathilda's answers get (`crosscheck.py`): zero verifies it, nonzero makes it wrong, and "
            "anything else leaves it unverified. Timeout 30 s per run.", "",
            "- **elementary** / **special**: a verified closed form, elementary or with special functions (Ei, li, "
            "erf, PolyLog, elliptic, hypergeometric, ...), classified by the same rule for every system. A special-function "
            "answer is counted, separately: on an integrand that is elementary it is correct but does not show the "
            "elementary antiderivative.",
            "- **decided**: closed, or proved non-elementary (Mathilda's certificate; FriCAS's unevaluated `integral`). "
            "A \"no\" on a case that is elementary by construction is a false \"no\" (Soundness below) and is not counted.",
            "- **Mathilda** is the Bronstein-family method of each class: RischTranscendental for R and T, Trager "
            "for A, Bronstein for M and B. **PMT** is ParallelMixedTower, the parallel Risch (Risch-Norman) "
            "heuristic over a mixed tower, run on every class; its \"no\" is a certificate too. **Automatic** is "
            "the full cascade.", ""]
    # the hand-written analysis kept beside the report (interpretation the tables
    # cannot carry: a check of the other system's failures, the timer floor, gaps)
    pre = os.path.splitext(a.out)[0] + ".preface.md"
    if os.path.exists(pre):
        head += open(pre).read().rstrip("\n").split("\n") + [""]
    S = soundness(allrows)
    sc = collections.Counter(k for k, _ in S)
    head.append("## Soundness")
    head.append("")
    if not S:
        head.append("No false \"no\" and no wrong answer from either system.")
    else:
        head.append("| finding | cases |")
        head.append("|---|---|")
        for k, v in sorted(sc.items()):
            head.append(f"| {k} | {v} |")
        head.append("")
        shown = collections.Counter()
        for k, r in sorted(S, key=lambda kr: ("fails" in kr[0], kr[0])):   # false "no" first
            shown[k] += 1
            if shown[k] <= 12:
                head.append(f"- {k}: `{r['id']}` `{r['f'][:140]}`")
            elif shown[k] == 13:
                head.append(f"- {k}: ... ({sc[k] - 12} more)")
    head.append("")
    if progress:
        head.append("## Mathilda against the previous saved run of each corpus")
        head.append("")
        head.append("| corpus | n | Mathilda closed | Mathilda decided | Automatic closed | Automatic decided |")
        head.append("|---|---|---|---|---|---|")
        for title, n, d0, d1, a0, a1 in progress:
            head.append(f"| {title} | {n} | {d0[0]} -> {d1[0]} | {d0[1]} -> {d1[1]} | {a0[0]} -> {a1[0]} | {a0[1]} -> {a1[1]} |")
        head.append("")
    oo = only_one(allrows)
    tail = ["## Decided by one system only (all corpora, by family)", "",
            "| family | Mathilda only | FriCAS only |", "|---|---|---|"]
    fams = sorted({f for _, f in oo})
    rows_oo = sorted(fams, key=lambda f: -(oo[("Mathilda only", f)] + oo[("FriCAS only", f)]))
    for f in rows_oo:
        if oo[("Mathilda only", f)] or oo[("FriCAS only", f)]:
            tail.append(f"| {f} | {oo[('Mathilda only', f)]} | {oo[('FriCAS only', f)]} |")
    tail.append("")
    mo = method_overlap(allrows)
    tail += ["## Within Mathilda: ParallelMixedTower against the Bronstein-family method (all corpora)", "",
             "Decided by one of the two only, by class and family.", "",
             "| class | family | PMT only | Bronstein-family only |", "|---|---|---|---|"]
    keys = sorted({(k[1], k[2]) for k in mo}, key=lambda cf: (-(mo[("PMT only",) + cf] + mo[("Bronstein-family only",) + cf]), cf))
    for cl, fam in keys:
        tail.append(f"| {cl} | {fam} | {mo[('PMT only', cl, fam)]} | {mo[('Bronstein-family only', cl, fam)]} |")
    tail.append("")
    fe = fricas_errors(allrows)
    tail += ["## Why FriCAS did not decide (all corpora)", "", "| reason | cases |", "|---|---|"]
    for k, v in fe.most_common(12):
        tail.append(f"| {k} | {v} |")
    tail.append("")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    open(a.out, "w").write("\n".join(head + L + tail) + "\n")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
