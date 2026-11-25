#!/usr/bin/env python3
"""latex_report.py -- a LaTeX/PDF report of one integration-shootout run.

    python3 tools/integration_shootout/latex_report.py --bench benchmarks/261-compact-integrals-x4 [--date D] [--no-pdf]

Reads the benchmark's corpus.tsv and results/<date>.json (the newest when no
--date) and writes <bench>/report/report.tex, compiled to report.pdf with
latexmk.  The report has
  * how every answer was verified, per system;
  * the tallies of every run -- each Mathilda method, FriCAS and the other
    systems -- with elementary and special-function answers counted apart, by the
    same rule for every system (translate.special_q on the system's own answer);
  * the head-to-head against FriCAS (elementary answers; special answers apart);
  * every hard failure and every answer verified only exactly;
  * an appendix with every integral of the corpus and every answer returned for
    it, with its status, time and verification.
Every formula is shown twice: verbatim, as the system returned it (so nothing
is lost in translation), and typeset.  The verbatim text is the full answer from
results/cache.jsonl (results/<date>.json keeps only its first 300 characters).
TeX comes from latex_tex.m: the integrand and Mathilda's and Wolfram's answers
are read as they are; FriCAS, Axiom, Maxima and sympy answers are first
translated into Mathilda syntax by crosscheck.to_wl, the translator the
cross-check uses.  A formula that cannot be typeset says why.
"""
import argparse
import collections
import json
import os
import re
import statistics
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import adapters as A  # noqa: E402
from tex_break import display, integral  # noqa: E402
from easy_classify import easy_ids  # noqa: E402
import crosscheck as X  # noqa: E402

SYSTEM = {"FriCAS": "fricas", "Axiom": "axiom", "Maxima": "maxima", "sympy": "sympy", "Wolfram": "wolfram"}

MATHILDA_METHODS = ["Automatic", "ParallelMixedTower", "RischTranscendental", "Trager", "Bronstein"]
EXTERNAL = ["FriCAS", "Axiom", "Maxima", "sympy", "Wolfram"]
DEDICATED = {"R": "RischTranscendental", "T": "RischTranscendental", "A": "Trager", "M": "Bronstein", "B": "Bronstein"}
CLASSNAME = {"R": "rational", "T": "transcendental (Log/Exp tower)", "A": "algebraic (one radical)",
             "M": "a radical over a Log/Exp tower", "B": "a radical below a Log/Exp kernel"}
STATUSES = ["ELEM", "SPECIAL", "CERT", "UNEVAL", "TIMEOUT", "ERROR", "UNVERIFIED", "WRONG"]
CLOSED = ("ELEM", "SPECIAL")


def tex_escape(s):
    rep = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_",
           "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    # a message name (Integrate::nonelem) may break after its "::"
    return "".join(rep.get(ch, ch) for ch in s).replace("::", r"::\allowbreak{}")


def load_corpus(bench):
    cases = {}
    for ln in open(os.path.join(bench, "corpus.tsv")).read().splitlines():
        q = ln.split("\t")
        if len(q) > 7 and q[0] == "CASE":
            cases[q[1]] = {"id": q[1], "class": q[2], "family": q[3], "expect": q[4], "f": q[5],
                           "F": None if q[6] == "None" else q[6], "points": q[7]}
    return cases


def load_rows(bench, date):
    res = os.path.join(bench, "results")
    if date is None:
        dates = sorted(f[:-5] for f in os.listdir(res) if re.match(r"\d{4}-\d\d-\d\d\.json$", f))
        date = dates[-1]
    return date, json.load(open(os.path.join(res, f"{date}.json")))["rows"]


def full_texts(bench, rows):
    """Replace each answer's 300-character prefix (results/<date>.json) by its full text from
    results/cache.jsonl: the last cached record of that system, method and case that extends
    it.  Marks raw["cut"] when the cached text itself stopped at the recording limit."""
    path = os.path.join(bench, "results", "cache.jsonl")
    texts = collections.defaultdict(list)
    for ln in (open(path) if os.path.exists(path) else ()):
        try:
            rec = json.loads(ln)
        except ValueError:
            continue
        k = (rec.get("ckey") or "").split("|")
        if len(k) >= 3 and rec.get("text"):
            texts[(k[0], k[1], k[2])].append(rec["text"])
    for r in rows:
        for lab, raw in (r.get("raw") or {}).items():
            pre = raw.get("text") or ""
            if not pre:
                continue
            key = ("mathilda", lab.split(":", 1)[1], r["id"]) if lab.startswith("Mathilda:") else \
                  (SYSTEM.get(lab, lab.lower()), "", r["id"])
            full = [t for t in texts.get(key, ()) if t.startswith(pre)]
            if full:
                raw["text"] = full[-1]
            # 4000: adapters.TRUNC's default, the recording limit before SHOOTOUT_TRUNC was raised
            raw["cut"] = len(raw["text"]) == 4000 or (not full and len(pre) == 300)


VERB_WIDTH = 112    # characters per verbatim line (\scriptsize \ttfamily on the text block)


def wrap_verbatim(text, width=VERB_WIDTH):
    """Hard-wrap an answer's text into lines of at most `width` characters, breaking
    after a space or an operator where one is near the end of the line.  listings
    can wrap a long line itself, but it builds the whole line in memory first, and
    a 30000-character answer on one line exhausts TeX (262, 263)."""
    out = []
    for line in text.split("\n"):
        while len(line) > width:
            cut = max(line.rfind(ch, width // 2, width) for ch in " +-*/,")
            cut = cut + 1 if cut > 0 else width
            out.append(line[:cut].rstrip())
            line = "  " + line[cut:].lstrip()
        out.append(line)
    return "\n".join(out)


def overview_page(rows, texs):
    """One page holding every integral of the corpus, typeset and scaled to fit.

    The integrals are set in a ragged paragraph box.  Its width is found by
    LaTeX itself: the box is set at a trial width, measured, and reset at the
    width that gives it the page's proportions (the area of a paragraph is
    nearly independent of its width), twice.  The box is then scaled by the
    smaller of textwidth/width and textheight/height, so it always fits the page
    and the formulas never overlap: they are only ever shrunk together."""
    order = sorted(rows, key=lambda r: ("RTAMB".index(r["class"]), r["id"]))
    from tex_break import split_terms, tokens

    def one(t):     # a sum is bracketed, as in the appendix's integral()
        return t if len(split_terms(tokens(t))) == 1 else r"\left(" + t + r"\right)"
    items = [r"\mbox{$\int " + one(texs[f"f|{r['id']}"]) + r"\,dx$}" for r in order if texs.get(f"f|{r['id']}")]
    if not items:
        return []
    setbox = (r"\sbox{\allintegrals}{\parbox{\allwd}{\raggedright\allintegralitems}}")
    # the page's proportions, and a width that gives the measured box those proportions
    refit = (r"\setlength{\allwd}{\fpeval{max(150, sqrt((\the\wd\allintegrals)*"
             r"((\the\ht\allintegrals)+(\the\dp\allintegrals))*(\the\textwidth)/(0.93*(\the\textheight))))}pt}")
    return [r"\begin{center}\textbf{Every integral of the corpus} (%d)\end{center}" % len(items),
            r"\newsavebox{\allintegrals}\newlength{\allwd}\setlength{\allwd}{1600pt}",  # a trial at 400pt is taller than TeX's 16383pt
            r"\newcommand{\allintegralitems}{%",
            "\n\\hspace{1.2em plus 1em} ".join(items) + "}",
            setbox, refit, setbox, refit, setbox,
            # never enlarged past 1.5 (a small corpus), never wider or taller than the page
            r"\begin{center}\scalebox{\fpeval{min(1.5, (\the\textwidth)/(\the\wd\allintegrals),"
            r"0.93*(\the\textheight)/((\the\ht\allintegrals)+(\the\dp\allintegrals)))}}"
            r"{\usebox{\allintegrals}}\end{center}",
            r"\clearpage",
            r"\sbox{\allintegrals}{}"]     # free it: a 1850-formula box held TeX's memory to the end (263)


def labels_present(rows):
    seen = set()
    for r in rows:
        seen.update(r["runs"])
    out = [f"Mathilda:{m}" for m in MATHILDA_METHODS if f"Mathilda:{m}" in seen]
    return out + [e for e in EXTERNAL if e in seen]


def short(lab):
    return lab.replace("Mathilda:", "").replace("ParallelMixedTower", "PMT").replace("RischTranscendental", "RischTr.")


def tex_batch(items, jobs=8, chunk=150):
    """{key: mathilda-syntax text} -> {key: TeX} through latex_tex.m, in parallel chunks.
    latex_tex.m writes a line for every key it reaches, so when a process dies on
    one expression (a crash, or the per-chunk time limit) that key is dropped and
    the chunk resumes after it: one bad expression costs only itself."""
    todo = [(k, v) for k, v in items.items() if "\t" not in v and "\n" not in v]
    if not todo:
        return {}
    d = tempfile.mkdtemp(prefix="latex-tex-")

    def run(part):
        out, n = {}, 0
        while part:
            n += 1
            fin, fout = os.path.join(d, f"in{id(part)}-{n}.tsv"), os.path.join(d, f"out{id(part)}-{n}.tsv")
            with open(fin, "w") as fh:
                fh.writelines(f"{k}\t{v}\n" for k, v in part)
            try:
                subprocess.run([A.MATHILDA, "-file", os.path.join(HERE, "latex_tex.m")],
                               env=dict(os.environ, IN=fin, OUT=fout), cwd=A.REPO, capture_output=True,
                               timeout=60 + 45 * len(part))
            except subprocess.TimeoutExpired:
                pass
            done = 0
            if os.path.exists(fout):
                for ln in open(fout).read().splitlines():
                    k, _, t = ln.partition("\t")
                    done += 1
                    if t:
                        out[k] = t
            part = part[done + 1:] if done < len(part) else []    # skip the key it died on
        return out

    res = {}
    with ThreadPoolExecutor(jobs) as ex:
        for got in ex.map(run, [todo[i:i + chunk] for i in range(0, len(todo), chunk)]):
            res.update(got)
    return res


def verification_text(lab, raw):
    nums = raw.get("nums") or []
    if isinstance(nums, str):
        nums = re.findall(r"ok|bad|na", nums)
    st = raw.get("status")
    det = raw.get("detail") or ""
    if st == "UNEVAL" and str(raw.get("cert")) == "True":
        return "certificate of non-elementarity (Integrate::nonelem)" + (f": {det}" if det else "")
    if st != "ANSWER":
        return {"UNEVAL": "returned unevaluated", "TIMEOUT": "no answer within the time limit",
                "ERROR": "error" + (f": {det}" if det else "")}.get(st, st)
    pts = ", ".join(nums) if nums else "none"
    if lab.startswith("Mathilda:"):
        own = "in Mathilda"
    elif lab in ("FriCAS", "Axiom"):
        own = f"complexNumeric in {lab}"
    else:
        own = f"in {lab}"
    v = raw.get("pzq", raw.get("xcheck"))
    verdict = {
        "zero": "PossibleZeroQ[D[F, x] - f] is True: verified",
        "nonzero": "PossibleZeroQ[D[F, x] - f] is False: wrong",
        "untranslatable": "the answer does not translate into Mathilda, so PossibleZeroQ cannot judge it: unverified",
        "cut": "the answer text was cut when the run was recorded, so PossibleZeroQ cannot judge it: unverified",
        None: "not judged by PossibleZeroQ"}.get(
        v, "PossibleZeroQ[D[F, x] - f] did not decide: unverified")
    return f"{verdict} (own three-point check, {own}: [{pts}])" + (f"; {det}" if det else "")


# The title and the opening paragraph of each benchmark's report.
_SEARCH = (r"They were found by a randomised search over class grammars (tools/integration\_shootout/gen\_compact.m, "
           r"plan {plan}): an integrand qualifies with Mathilda's verified answer as reference, or, when Mathilda neither "
           r"answers, answers wrongly, nor answers elementarily, with a verified answer of FriCAS, Axiom or Maxima "
           r"(families \texttt{{-fricas}}, \texttt{{-axiom}}, \texttt{{-maxima}}); Wolfram is never consulted. ")
_HARD = (r"Every integrand is outside the \emph{{easy}} classes of easy.m: it is not a rational function, no "
         r"linear-fractional, Euler or chained substitution makes it rational, and Integrate with "
         r"\texttt{{Method -> \char`\"DerivativeDivides\char`\"}} does not close it. ")
ABOUT = {
    "256-integration-shootout": {
        "title": "the integration shootout",
        "intro": lambda n: (f"The corpus has {n} indefinite integrals $\\int f\\,dx$: the shootout's "
                            r"families, generated by gen.m (seeds 1--3): rational functions with repeated poles and RootSum logarithms, "
                            r"Log/Exp towers, Trager's algebraic stress corpus, Bronstein's mixed families and radicals "
                            r"below a logarithm or an exponential, each with perturbed variants; elementary, provably "
                            r"non-elementary, or unknown.")},
    "257-compact-integrals": {
        "title": "compact integrands, integrals at least twice as large",
        "intro": lambda n: (f"The corpus has {n} indefinite integrals $\\int f\\,dx$. Each integrand $f$ is compact, and has "
                            r"an elementary antiderivative $F$ whose compact form is at least twice as large: "
                            r"$2\,\mathrm{LeafCount}(f) \le \mathrm{LeafCount}(F)$. " + _SEARCH.format(plan=257))},
    "259-transcendental-decision": {
        "title": "the transcendental decision corpus",
        "intro": lambda n: (f"The corpus has {n} indefinite integrals $\\int f\\,dx$, all transcendental elementary "
                            r"functions with no radical of $x$: the families of the completeness campaign's hunts. A "
                            r"family ending in \texttt{-elem} is $f = F'$, elementary by construction; one ending in "
                            r"\texttt{-rand} is random, elementary or not.")},
    "260-radical-bottom": {
        "title": "a radical at the bottom",
        "intro": lambda n: (f"The corpus has {n} indefinite integrals $\\int f\\,dx$, all with a simple radical of $x$ "
                            r"below a logarithm or an exponential. A family ending in \texttt{-elem} is $f = F'$, "
                            r"elementary by construction; the others are random or perturbed, elementary or not.")},
    "261-compact-integrals-x4": {
        "title": "compact integrands, integrals at least four times larger",
        "intro": lambda n: (f"The corpus has {n} indefinite integrals $\\int f\\,dx$. Each integrand $f$ is compact, and has "
                            r"an elementary antiderivative $F$ whose compact form is more than four times as large: "
                            r"$4\,\mathrm{LeafCount}(f) < \mathrm{LeafCount}(F)$. " + _SEARCH.format(plan=261) +
                            r"The integrands of benchmarks 256, 257, 259 and 260 are excluded.")},
    "262-hard-integrals": {
        "title": "integrals outside the easy classes",
        "intro": lambda n: (f"The corpus has {n} indefinite integrals $\\int f\\,dx$, each with an elementary antiderivative. " +
                            _HARD.format() + r"Classes T, A and B need $\mathrm{LeafCount}(F) \ge 2\,\mathrm{LeafCount}(f)$; "
                            r"M, whose answers are built on a unit $\log(P + \sqrt{P^2 - c e^h})$, only "
                            r"$\mathrm{LeafCount}(F) > \mathrm{LeafCount}(f)$. " + _SEARCH.format(plan=262) +
                            r"The integrands of benchmarks 256, 257 and 259--261 are excluded.")},
    "264-mixed-integrals": {
        "title": "a radical on top of a Log/Exp tower",
        "intro": lambda n: (f"The corpus has {n} indefinite integrals $\\int f\\,dx$, all of class M: a radical whose "
                            r"radicand holds an exponential of $x$, each with an elementary antiderivative. " +
                            _HARD.format() + r"They are built on a unit of the tower: with $g = P^2 - c\,e^h$, "
                            r"$(P' - h' P/2)/\sqrt{g}$ integrates to $\log(P + \sqrt{g}) - h/2$, and each integrand "
                            r"is drawn from that shape and its perturbations (grammar M4), with "
                            r"$\mathrm{LeafCount}(F) > \mathrm{LeafCount}(f)$. " + _SEARCH.format(plan=264) +
                            r"The integrands of benchmarks 256, 257 and 259--263 are excluded.")},
    "263-hard-integrals-ii": {
        "title": "more integrals outside the easy classes",
        "intro": lambda n: (f"The corpus has {n} indefinite integrals $\\int f\\,dx$, each with an elementary antiderivative. " +
                            _HARD.format() + r"Classes T and A need $\mathrm{LeafCount}(F) \ge 2\,\mathrm{LeafCount}(f)$, "
                            r"M only $\mathrm{LeafCount}(F) > \mathrm{LeafCount}(f)$; B has no size condition (its B4 "
                            r"integrands are derivatives of products, so $F$ is the smaller). " + _SEARCH.format(plan=263) +
                            r"The integrands of benchmarks 256, 257 and 259--262 are excluded.")},
}


def about_of(name):
    return ABOUT.get(name, {"title": "integration shootout",
                            "intro": lambda n: f"The corpus has {n} indefinite integrals $\\int f\\,dx$."})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", required=True)
    ap.add_argument("--date", default=None)
    ap.add_argument("--no-pdf", action="store_true")
    ap.add_argument("--limit", type=int, default=0, help="first N cases per class only (a quick check of the document)")
    ap.add_argument("--out", default=None, help="report folder (default <bench>/report)")
    ap.add_argument("--include-easy", action="store_true",
                    help="keep the cases easy.tsv classes as easy (left out by default)")
    a = ap.parse_args()
    bench = os.path.abspath(a.bench)
    cases = load_corpus(bench)
    date, rows = load_rows(bench, a.date)
    easy = set() if a.include_easy else easy_ids(bench)
    n_easy = sum(1 for r in rows if r["id"] in cases and r["id"] in easy)
    rows = [r for r in rows if r["id"] in cases and r["id"] not in easy]
    full_texts(bench, rows)
    if a.limit:
        per = collections.Counter()
        rows = [r for r in rows if (per.update([r["class"]]) or per[r["class"]] <= a.limit)]
    labs = labels_present(rows)
    name = os.path.basename(bench)
    about = about_of(name)

    # ---------------- tallies ----------------
    tally = {lab: collections.Counter() for lab in labs}
    times = {lab: [] for lab in labs}
    bycls = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
    for r in rows:
        for lab, cls in r["runs"].items():
            if lab not in tally:
                continue
            tally[lab][cls] += 1
            tally[lab]["n"] += 1
            bycls[r["class"]][lab][cls] += 1
            bycls[r["class"]][lab]["n"] += 1
            t = (r.get("raw") or {}).get(lab, {}).get("time")
            if cls in CLOSED and t not in (None, "", "None"):
                times[lab].append(float(t))

    def best_mathilda(r):
        cl = [r["runs"].get(f"Mathilda:{m}") for m in MATHILDA_METHODS]
        return "ELEM" if "ELEM" in cl else ("SPECIAL" if "SPECIAL" in cl else ("CERT" if "CERT" in cl else "OTHER"))

    def h2h(get_m):
        c = collections.Counter()
        for r in rows:
            if "FriCAS" not in r["runs"]:
                continue
            m, f = get_m(r), r["runs"]["FriCAS"]
            me, fe = m == "ELEM", f == "ELEM"
            c[("both" if me and fe else "Mathilda only" if me else "FriCAS only" if fe else "neither")] += 1
            if m == "SPECIAL":
                c["Mathilda special"] += 1
            if f == "SPECIAL":
                c["FriCAS special"] += 1
            c["n"] += 1
        return c

    h_best = h2h(best_mathilda)
    h_ded = h2h(lambda r: r["runs"].get(f"Mathilda:{DEDICATED[r['class']]}"))
    hard = []
    for r in rows:
        for v in r.get("verdicts", []):
            if v.split()[0] in ("WRONG", "FALSE_CERT", "DISAGREE"):
                hard.append((r["id"], v))
        for lab in labs:
            if not lab.startswith("Mathilda:") and r["runs"].get(lab) == "WRONG":
                hard.append((r["id"], f"{lab} answer rejected by PossibleZeroQ[D[F, x] - f] in Mathilda"))
    xtab = collections.defaultdict(collections.Counter)
    for r in rows:
        for lab in labs:
            raw = (r.get("raw") or {}).get(lab) or {}
            v = raw.get("pzq", raw.get("xcheck"))
            if v is not None and raw.get("status") == "ANSWER":
                nums = raw.get("nums") or []
                nums = nums if isinstance(nums, list) else re.findall(r"ok|bad|na", nums)
                own = "ok" if nums and all(x == "ok" for x in nums) else ("bad" if "bad" in nums else "na")
                v = v if v in ("zero", "nonzero", "untranslatable", "cut") else "undecided"
                xtab[(lab, own)][v] += 1
    exact_only = [(r["id"], (r.get("raw") or {}).get("FriCAS", {}).get("detail", "")) for r in rows
                  if "point check failed" in ((r.get("raw") or {}).get("FriCAS", {}).get("detail") or "")]

    # ---------------- TeX of integrands and Mathilda answers ----------------
    items, untyped = {}, {}
    for r in rows:
        items[f"f|{r['id']}"] = cases[r["id"]]["f"]
        for lab in labs:
            raw = (r.get("raw") or {}).get(lab)
            if not (raw and raw.get("status") == "ANSWER" and raw.get("text")) or raw.get("cut"):
                continue
            if lab.startswith("Mathilda:") or lab == "Wolfram":
                items[f"{lab}|{r['id']}"] = raw["text"]
            else:
                try:
                    items[f"{lab}|{r['id']}"] = X.to_wl(SYSTEM[lab], raw["text"])
                except Exception as e:  # noqa: BLE001 -- any parse failure: not typeset, and said so
                    untyped[f"{lab}|{r['id']}"] = "the answer does not translate into Mathilda syntax" + (
                        f" ({e})" if isinstance(e, X.Untranslatable) and str(e) else "")
    texs = tex_batch(items)

    # ---------------- the document ----------------
    L = [r"\documentclass[10pt]{article}",
         r"\usepackage[a4paper,margin=18mm]{geometry}",
         r"\usepackage{amsmath,amssymb,booktabs,longtable,listings,xcolor,graphicx}",
         r"\usepackage[hidelinks]{hyperref}",
         r"\usepackage{xfp}",
         r"\lstset{basicstyle=\ttfamily\scriptsize,breaklines=true,breakatwhitespace=false,columns=fullflexible,"
         r"keepspaces=true,showstringspaces=false,frame=none,xleftmargin=1em,"
         # a break point after every digit, so a 150-digit integer can wrap
         r"literate=" + " ".join("{%d}{{%d}\\allowbreak}1" % (d, d) for d in range(10)) + "}",
         r"\setlength{\parindent}{0pt}\setlength{\parskip}{3pt}\setlength{\emergencystretch}{3em}",
         r"\allowdisplaybreaks",
         r"\title{Benchmark " + tex_escape(name) + ": " + about["title"] + "}",
         r"\author{Mathilda integration shootout}",
         r"\date{Run " + date + "}",
         r"\begin{document}", r"\maketitle", r"\tableofcontents", r"\clearpage"]
    L += overview_page(rows, texs)

    n = len(rows)
    L += [r"\section{The benchmark}",
          about["intro"](n),
          (f"{n_easy} further cases of the corpus are \\emph{{easy}} and are left out of this report: "
           r"easy.m (tools/integration\_shootout) classes each of them as a rational function, an integrand "
           r"that a linear-fractional, Euler or chained substitution makes rational, or one the "
           r"DerivativeDivides method closes (easy.tsv lists the class of every case)." if n_easy else ""), "",
          r"\begin{center}\begin{tabular}{lll}\toprule class & integrands & cases\\\midrule"]
    cc = collections.Counter(r["class"] for r in rows)
    for k in "RTAMB":
        if cc[k]:
            L.append(f"{k} & {tex_escape(CLASSNAME[k])} & {cc[k]}\\\\")
    L += [r"\bottomrule\end{tabular}\end{center}"]

    L += [r"\section{How each answer was verified}",
          r"Every answer of every system is judged by one test, in Mathilda: "
          r"\texttt{PossibleZeroQ[D[F, x] - f]} (60~s), the test the corpus itself was built under "
          r"(Section~\ref{sec:xcheck}). \textbf{True}: verified. \textbf{False}: \textbf{wrong}. Anything else "
          r"is \textbf{unverified}, never verified on numeric evidence: the test did not decide, or the answer "
          r"does not translate into Mathilda, or its text was cut when the run was recorded. A True counts only "
          r"where $D[F,x]-f$ evaluates to a number at two test points, since \texttt{PossibleZeroQ} answers True "
          r"for an expression it cannot evaluate.",
          r"\begin{itemize}",
          r"\item \textbf{Mathilda} (every method) and \textbf{Wolfram}: the answer is read as it is. "
          r"A \emph{certificate} is Mathilda's proof that no elementary antiderivative exists "
          r"(\texttt{Integrate::nonelem}); on an integrand that is elementary by construction it is a false "
          r"``no'' and a hard failure.",
          r"\item \textbf{FriCAS, Axiom} (Tim Daly's distribution, the system FriCAS forked), \textbf{Maxima, "
          r"sympy}: the answer is translated from the system's output into Mathilda syntax first. FriCAS's "
          r"\texttt{rootOf(p, v)}, valid for every root of $p$, becomes \texttt{Root[p \&, 1]}. The translator was "
          r"checked both ways: answers multiplied by $1 + x/10^6$ all come out nonzero.",
          r"\item \textbf{Each system's own check} is still run and shown after the verdict, as a secondary "
          r"note: $D[F,x]-f$ at three rational points to 40 digits (FriCAS and Axiom: \texttt{complexNumeric}, "
          r"and \texttt{normalize(D(r, x) - f)} for answers it cannot evaluate). It decides nothing; where it "
          r"disagrees with \texttt{PossibleZeroQ}, the verdict is \texttt{PossibleZeroQ}'s.",
          r"\item A verified answer is \textbf{elementary} or \textbf{special} by the functions it contains "
          r"(Ei, li, erf, Si, Ci, polylogarithms, elliptic, hypergeometric, Bessel, Gamma, \ldots), by the same "
          r"list for every system. A special-function answer is correct and is counted, but separately: it does not "
          r"exhibit the elementary antiderivative these integrands have.",
          r"\item An \textbf{error} is shown with the message the system printed.",
          r"\end{itemize}",
          r"Time limit: 30 s per run."]

    L += [r"\section{Results}", r"\subsection{Every run}",
          r"\begin{center}\small\begin{tabular}{l" + "r" * (len(STATUSES) + 2) + r"}\toprule",
          "run & n & " + " & ".join(s.lower() for s in STATUSES) + r" & median s\\\midrule"]
    for lab in labs:
        c = tally[lab]
        med = f"{statistics.median(times[lab]):.2f}" if times[lab] else "--"
        L.append(tex_escape(short(lab) if lab.startswith("Mathilda:") else lab) +
                 (" (Mathilda)" if lab.startswith("Mathilda:") else "") +
                 f" & {c['n']} & " + " & ".join(str(c[s]) for s in STATUSES) + f" & {med}\\\\")
    L += [r"\bottomrule\end{tabular}\end{center}",
          r"Each Mathilda method runs on the classes it applies to: Automatic and ParallelMixedTower (PMT) on all, "
          r"RischTranscendental on R and T, Trager on A, Bronstein on M and B."]

    L += [r"\subsection{By class: elementary / special / decided}"]
    for k in "RTAMB":
        if not cc[k]:
            continue
        L += [f"\\paragraph{{Class {k}: {tex_escape(CLASSNAME[k])} ({cc[k]} cases)}}",
              r"\begin{center}\small\begin{tabular}{lrrrr}\toprule run & elementary & special & certificate & "
              r"other\\\midrule"]
        for lab in labs:
            c = bycls[k][lab]
            if not c["n"]:
                continue
            oth = c["n"] - c["ELEM"] - c["SPECIAL"] - c["CERT"]
            L.append(tex_escape(short(lab)) + f" & {c['ELEM']} & {c['SPECIAL']} & {c['CERT']} & {oth}\\\\")
        L += [r"\bottomrule\end{tabular}\end{center}"]

    def h2h_table(h, title):
        return [f"\\paragraph{{{title}}}",
                r"\begin{center}\begin{tabular}{lr}\toprule & cases\\\midrule",
                f"both elementary & {h['both']}\\\\", f"Mathilda only & {h['Mathilda only']}\\\\",
                f"FriCAS only & {h['FriCAS only']}\\\\", f"neither & {h['neither']}\\\\\\midrule",
                f"Mathilda special-function answers & {h['Mathilda special']}\\\\",
                f"FriCAS special-function answers & {h['FriCAS special']}\\\\",
                r"\bottomrule\end{tabular}\end{center}"]

    L += [r"\subsection{Against FriCAS}",
          r"Elementary answers, case by case. ``Mathilda best'' is any Mathilda method; ``dedicated'' is the "
          r"decision procedure of the class (RischTranscendental, Trager or Bronstein)."]
    L += h2h_table(h_best, "Mathilda best against FriCAS")
    L += h2h_table(h_ded, "Mathilda dedicated method against FriCAS")

    L += [r"\subsection{Hard failures}"]
    if hard:
        L += [r"\begin{longtable}{ll}\toprule case & verdict\\\midrule"]
        L += [f"\\hyperlink{{case-{cid}}}{{{tex_escape(cid)}}} & {tex_escape(v)}\\\\" for cid, v in hard]
        L += [r"\bottomrule\end{longtable}"]
    else:
        L.append("None.")
    L += [r"\subsection{PossibleZeroQ verdicts}\label{sec:xcheck}",
          r"Every answer, judged by \texttt{PossibleZeroQ[D[F, x] - f]} in Mathilda, by the result of the "
          r"system's own three-point check (ok: passed at every point; bad: failed at a point; na: no point "
          r"could be evaluated). A row where the own check is ok and the verdict is not zero, or the other way "
          r"round, is where the two disagree."]
    if xtab:
        L += [r"\begin{center}\small\begin{tabular}{llrrrrrr}\toprule run & own check & n & zero & nonzero & "
              r"undecided & untranslatable & cut\\\midrule"]
        for lab in labs:
            for own in ("ok", "bad", "na"):
                c = xtab.get((lab, own))
                if c:
                    L.append(f"{tex_escape(short(lab))} & {own} & {sum(c.values())} & {c['zero']} & {c['nonzero']} & "
                             f"{c['undecided']} & {c['untranslatable']} & {c['cut']}\\\\")
        L += [r"\bottomrule\end{tabular}\end{center}"]
    else:
        L.append("None.")
    L += [r"\subsection{FriCAS answers verified only exactly}\label{sec:exact}"]
    if exact_only:
        L += [r"\begin{longtable}{ll}\toprule case & detail\\\midrule"]
        L += [f"\\hyperlink{{case-{cid}}}{{{tex_escape(cid)}}} & {tex_escape(d)}\\\\" for cid, d in exact_only]
        L += [r"\bottomrule\end{longtable}"]
    else:
        L.append("None.")

    # ---------------- appendix ----------------
    L += [r"\clearpage", r"\appendix", r"\section{Every integral and every answer}"]
    order = sorted(rows, key=lambda r: ("RTAMB".index(r["class"]), r["id"]))
    for r in order:
        c = cases[r["id"]]
        L.append(f"\\subsubsection*{{\\hypertarget{{case-{r['id']}}}{{{tex_escape(r['id'])}}} "
                 f"\\normalfont\\small(class {c['class']}, family {tex_escape(c['family'])})}}")
        ft = texs.get(f"f|{r['id']}")
        L += [r"\begin{lstlisting}", wrap_verbatim("Integrate[" + c["f"] + ", x]"), r"\end{lstlisting}"]
        L.append(integral(ft) if ft else r"\par\emph{(not typeset: TeXForm failed)}")
        L.append("")
        for lab in labs:
            raw = (r.get("raw") or {}).get(lab)
            cls = r["runs"].get(lab)
            if raw is None or cls is None:
                continue
            t = raw.get("time")
            ts = f", {float(t):.2f} s" if t not in (None, "", "None") else ""
            st = raw.get("status")
            if st == "ERROR":       # the system's own error message, as it printed it
                what = r"\emph{Error message:} " + tex_escape(raw.get("detail") or "none recorded")
            elif st in ("UNEVAL", "TIMEOUT") and not (st == "UNEVAL" and str(raw.get("cert")) == "True"):
                what = tex_escape(verification_text(lab, raw))
            else:
                what = r"\emph{Verdict:} " + tex_escape(verification_text(lab, raw))
            L.append(f"\\textbf{{{tex_escape(short(lab) if lab.startswith('Mathilda:') else lab)}}}"
                     + (" (Mathilda)" if lab.startswith("Mathilda:") else "")
                     + f": {cls}{ts}. {what}")
            L.append("")    # every run is its own paragraph
            if raw.get("status") == "ANSWER" and raw.get("text"):
                k = f"{lab}|{r['id']}"
                L += [r"\begin{lstlisting}", wrap_verbatim(raw["text"]), r"\end{lstlisting}"]
                if raw.get("cut"):
                    L.append(r"\emph{(the answer text was cut at %d characters when the run was recorded, "
                             r"so it is not typeset)}" % len(raw["text"]))
                elif texs.get(k):
                    L.append(display("F", texs[k]))
                else:
                    L.append(r"\emph{(not typeset: %s)}" % tex_escape(untyped.get(k, "TeXForm failed")))
                L.append("")
    L += [r"\end{document}"]

    outdir = a.out or os.path.join(bench, "report")
    os.makedirs(outdir, exist_ok=True)
    tex = os.path.join(outdir, "report.tex")
    with open(tex, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"wrote {tex}: {n} cases, {len(labs)} runs, {len(texs)} typeset expressions")
    if not a.no_pdf:
        p = subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "report.tex"],
                           cwd=outdir, capture_output=True, text=True)
        log = os.path.join(outdir, "report.log")
        if p.returncode != 0 and os.path.exists(log) and "TeX capacity exceeded" in open(log, errors="replace").read():
            # pdfTeX's main memory is fixed (5M words); LuaTeX allocates as it goes.
            # 262 and 263, with every answer typeset, need the latter.
            print("pdfTeX ran out of memory; compiling with LuaLaTeX", flush=True)
            subprocess.run(["latexmk", "-C"], cwd=outdir, capture_output=True)
            p = subprocess.run(["latexmk", "-lualatex", "-interaction=nonstopmode", "-halt-on-error", "report.tex"],
                               cwd=outdir, capture_output=True, text=True)
        pdf = os.path.join(outdir, "report.pdf")
        print("pdf:", pdf if p.returncode == 0 and os.path.exists(pdf) else "FAILED -- see report.log")
        if p.returncode != 0:
            sys.exit(1)


if __name__ == "__main__":
    main()
