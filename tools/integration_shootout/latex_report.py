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
Answers are shown as the system returned them (verbatim, so nothing is lost
in translation), and typeset as well when Mathilda's TeXForm of them is short.
TeX for the Mathilda answers and the integrands comes from latex_tex.m.
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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import adapters as A  # noqa: E402
from tex_break import display  # noqa: E402
from easy_classify import easy_ids  # noqa: E402

MATHILDA_METHODS = ["Automatic", "ParallelMixedTower", "RischTranscendental", "Trager", "Bronstein"]
EXTERNAL = ["FriCAS", "Axiom", "Maxima", "sympy", "Wolfram"]
DEDICATED = {"R": "RischTranscendental", "T": "RischTranscendental", "A": "Trager", "M": "Bronstein", "B": "Bronstein"}
CLASSNAME = {"R": "rational", "T": "transcendental (Log/Exp tower)", "A": "algebraic (one radical)",
             "M": "a radical over a Log/Exp tower", "B": "a radical below a Log/Exp kernel"}
STATUSES = ["ELEM", "SPECIAL", "CERT", "UNEVAL", "TIMEOUT", "ERROR", "UNVERIFIED", "WRONG"]
CLOSED = ("ELEM", "SPECIAL")
TYPESET_MAX = 3000         # longest TeX typeset in the appendix (the raw text is always shown);
                           # tex_break.display wraps a long one over multline* lines


def tex_escape(s):
    rep = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_",
           "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(rep.get(ch, ch) for ch in s)


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


def labels_present(rows):
    seen = set()
    for r in rows:
        seen.update(r["runs"])
    out = [f"Mathilda:{m}" for m in MATHILDA_METHODS if f"Mathilda:{m}" in seen]
    return out + [e for e in EXTERNAL if e in seen]


def short(lab):
    return lab.replace("Mathilda:", "").replace("ParallelMixedTower", "PMT").replace("RischTranscendental", "RischTr.")


def tex_batch(items):
    """{key: mathilda-syntax text} -> {key: TeX} through latex_tex.m"""
    if not items:
        return {}
    with tempfile.TemporaryDirectory() as d:
        fin, fout = os.path.join(d, "in.tsv"), os.path.join(d, "out.tsv")
        with open(fin, "w") as fh:
            for k, v in items.items():
                if "\t" not in v and "\n" not in v:
                    fh.write(f"{k}\t{v}\n")
        subprocess.run([A.MATHILDA, "-file", os.path.join(HERE, "latex_tex.m")],
                       env=dict(os.environ, IN=fin, OUT=fout), cwd=A.REPO, capture_output=True,
                       timeout=12 * 3600)
        out = {}
        if os.path.exists(fout):
            for ln in open(fout).read().splitlines():
                k, _, t = ln.partition("\t")
                out[k] = t
        return out


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
        how = "D[F, x] - f at the three rational points, 40 digits, relative 1e-20, in Mathilda"
    elif lab in ("FriCAS", "Axiom"):
        how = f"complexNumeric(eval(D(r, x) - f)) at the three points, 40 digits, in {lab}"
    else:
        how = "the derivative at the three points, in the system itself"
    xv = raw.get("xcheck")
    if xv is not None:
        det = (det + "; " if det else "") + {
            "zero": "own check failed; translated into Mathilda, PossibleZeroQ[D[F, x] - f] is True: verified",
            "nonzero": "own check failed; translated into Mathilda, PossibleZeroQ[D[F, x] - f] is False: wrong",
            "untranslatable": "own check failed; the answer does not translate into Mathilda: unverified"}.get(
            xv, "own check failed; PossibleZeroQ[D[F, x] - f] in Mathilda did not decide: unverified")
    return f"{how}: [{pts}]" + (f"; {det}" if det else "")


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
    if a.limit:
        per = collections.Counter()
        rows = [r for r in rows if (per.update([r["class"]]) or per[r["class"]] <= a.limit)]
    labs = labels_present(rows)
    name = os.path.basename(bench)

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
                hard.append((r["id"], f"{lab} answer fails its own check and PossibleZeroQ in Mathilda"))
    xtab = collections.defaultdict(collections.Counter)
    for r in rows:
        for lab in labs:
            raw = (r.get("raw") or {}).get(lab) or {}
            if "xcheck" in raw:
                nums = raw.get("nums") or []
                own = "bad" if "bad" in (nums if isinstance(nums, list) else re.findall(r"ok|bad|na", nums)) else "na"
                v = raw["xcheck"] if raw["xcheck"] in ("zero", "nonzero", "untranslatable") else "undecided"
                xtab[(lab, own)][v] += 1
    exact_only = [(r["id"], (r.get("raw") or {}).get("FriCAS", {}).get("detail", "")) for r in rows
                  if "point check failed" in ((r.get("raw") or {}).get("FriCAS", {}).get("detail") or "")]

    # ---------------- TeX of integrands and Mathilda answers ----------------
    items = {}
    for r in rows:
        items[f"f|{r['id']}"] = cases[r["id"]]["f"]
        for lab in labs:
            raw = (r.get("raw") or {}).get(lab)
            if lab.startswith("Mathilda:") and raw and raw.get("status") == "ANSWER" and raw.get("text"):
                if len(raw["text"]) < 6000:
                    items[f"{lab}|{r['id']}"] = raw["text"]
    texs = tex_batch(items)

    # ---------------- the document ----------------
    L = [r"\documentclass[10pt]{article}",
         r"\usepackage[a4paper,margin=18mm]{geometry}",
         r"\usepackage{amsmath,amssymb,booktabs,longtable,listings,xcolor,graphicx}",
         r"\usepackage[hidelinks]{hyperref}",
         r"\lstset{basicstyle=\ttfamily\scriptsize,breaklines=true,breakatwhitespace=false,columns=fullflexible,"
         r"keepspaces=true,showstringspaces=false,frame=none,xleftmargin=1em}",
         r"\setlength{\parindent}{0pt}\setlength{\parskip}{3pt}",
         r"\allowdisplaybreaks",
         r"\title{Benchmark " + tex_escape(name) + r": compact integrands, integrals at least four times larger}",
         r"\author{Mathilda integration shootout}",
         r"\date{Run " + date + "}",
         r"\begin{document}", r"\maketitle", r"\tableofcontents", r"\clearpage"]

    n = len(rows)
    L += [r"\section{The benchmark}",
          f"The corpus has {n} indefinite integrals $\\int f\\,dx$. Each integrand $f$ is compact, and has an elementary "
          r"antiderivative $F$ whose compact form is more than four times as large: "
          r"$4\,\mathrm{LeafCount}(f) < \mathrm{LeafCount}(F)$. They were found by a randomised search over class "
          r"grammars (tools/integration\_shootout/gen\_compact.m, plan 261): an integrand qualifies with Mathilda's "
          r"verified answer as reference, or, when Mathilda neither answers, answers wrongly, nor answers elementarily, "
          r"with FriCAS's verified answer (families \texttt{-fricas}); Wolfram is never consulted. "
          r"The integrands of benchmarks 256, 257, 259 and 260 are excluded.",
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
          r"Every system first verifies its own answer. An answer that fails its own check, or that its own check "
          r"cannot evaluate, is then translated into Mathilda and judged by the test Mathilda's answers are generated "
          r"under, \texttt{PossibleZeroQ[D[F, x] - f]} (Section~\ref{sec:xcheck}). So \emph{wrong} means the same "
          r"for every system: the answer fails its own check \emph{and} \texttt{PossibleZeroQ}.",
          r"\begin{itemize}",
          r"\item \textbf{Mathilda} (every method): the derivative $D[F,x]-f$ is evaluated at three rational points "
          r"at 40 significant digits and must be below $10^{-20}$ relative to $f$ at every point. "
          r"A \emph{certificate} is Mathilda's proof that no elementary antiderivative exists "
          r"(\texttt{Integrate::nonelem}); on this corpus every integrand is elementary, so a certificate is a false "
          r"``no'' and a hard failure.",
          r"\item \textbf{FriCAS}: \texttt{complexNumeric(eval(D(r, x) - f, x, p), 40)} at the same three points, "
          r"in FriCAS. When a point cannot be evaluated, or when the point check fails but "
          r"\texttt{normalize(D(r, x) - f)} is exactly $0$, the answer is counted as verified (exactly); the second "
          r"case, a branch choice at the sample points, is listed in Section~\ref{sec:exact}.",
          r"\item \textbf{Axiom} (Tim Daly's distribution, the system FriCAS forked): the same check as FriCAS, run by "
          r"\texttt{AXIOMsys}.",
          r"\item \textbf{Maxima, sympy}: the derivative at the three points, in the system itself.",
          r"\item A verified answer is \textbf{elementary} or \textbf{special} by the functions it contains "
          r"(Ei, li, erf, Si, Ci, polylogarithms, elliptic, hypergeometric, Bessel, Gamma, \ldots), by the same "
          r"list for every system. A special-function answer is correct and is counted, but separately: it does not "
          r"exhibit the elementary antiderivative these integrands have.",
          r"\item \textbf{Cross-check in Mathilda}: an answer whose own check failed at a point (``bad'') or "
          r"evaluated at none (``na'') is parsed from the system's output and judged by "
          r"\texttt{PossibleZeroQ[D[F, x] - f]} (60~s). True: verified (a branch or an evaluation limit of the "
          r"system's own check). False: \textbf{wrong}. Undecided, or an answer that does not translate: "
          r"\textbf{unverified}, never wrong. The translator was checked on answers that passed their own checks: "
          r"every one of them is zero in Mathilda. FriCAS's \texttt{rootOf(p, v)}, valid for every root of $p$, "
          r"becomes \texttt{Root[p \&, 1]}.",
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
    L += [r"\subsection{Cross-check in Mathilda}\label{sec:xcheck}",
          r"Answers that failed (bad) or could not complete (na) their own system's check, judged by "
          r"\texttt{PossibleZeroQ[D[F, x] - f]} in Mathilda."]
    if xtab:
        L += [r"\begin{center}\small\begin{tabular}{llrrrrr}\toprule run & own check & n & zero & nonzero & "
              r"undecided & untranslatable\\\midrule"]
        for lab in labs:
            for own in ("bad", "na"):
                c = xtab.get((lab, own))
                if c:
                    L.append(f"{tex_escape(short(lab))} & {own} & {sum(c.values())} & {c['zero']} & {c['nonzero']} & "
                             f"{c['undecided']} & {c['untranslatable']}\\\\")
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
        if ft:
            L.append(display("f", ft))
        L += [r"\begin{lstlisting}", "f = " + c["f"], r"\end{lstlisting}"]
        for lab in labs:
            raw = (r.get("raw") or {}).get(lab)
            cls = r["runs"].get(lab)
            if raw is None or cls is None:
                continue
            t = raw.get("time")
            ts = f", {float(t):.2f} s" if t not in (None, "", "None") else ""
            L.append(f"\\textbf{{{tex_escape(short(lab) if lab.startswith('Mathilda:') else lab)}}}"
                     + (" (Mathilda)" if lab.startswith("Mathilda:") else "")
                     + f": {cls}{ts}. \\emph{{Verified:}} {tex_escape(verification_text(lab, raw))}")
            if raw.get("status") == "ANSWER" and raw.get("text"):
                at = texs.get(f"{lab}|{r['id']}")
                if at and len(at) <= TYPESET_MAX:
                    L.append(display("F", at))
                L += [r"\begin{lstlisting}", raw["text"], r"\end{lstlisting}"]
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
        pdf = os.path.join(outdir, "report.pdf")
        print("pdf:", pdf if p.returncode == 0 and os.path.exists(pdf) else "FAILED -- see report.log")
        if p.returncode != 0:
            sys.exit(1)


if __name__ == "__main__":
    main()
