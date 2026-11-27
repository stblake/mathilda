#!/usr/bin/env python3
"""compare_pdf.py -- a LaTeX/PDF listing of every integral of a Goursat stress
benchmark with the answer each requested Integrate method returned for it.

    python3 tools/goursat/compare_pdf.py --bench <dir> [--methods GeneralisedGoursat,Trager] [--out <dir>]

Reads <bench>/corpus.tsv, the newest <bench>/results/*.json (the per-case
status of every run, as goursat_stress.py judged it: PossibleZeroQ[D[F, x] - f])
and <bench>/results/cache.jsonl (the full answer text).  Run goursat_stress.py
first with the same --methods.  The document has

  * a summary: per family, the status counts of each method, and how often the
    two methods agree (both close, only one closes, neither);
  * every case, grouped by family: the integral, then each method's status,
    time, LeafCount and answer -- typeset (through latex_tex.m), with the raw
    text beside it when it cannot be typeset or is very long.
"""
import argparse
import collections
import glob
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "tools", "integration_shootout"))
sys.path.insert(0, HERE)
from tex_break import display, integral  # noqa: E402
from latex_report import tex_escape, tex_batch  # noqa: E402
import goursat_stress as GS  # noqa: E402

CLOSED = ("ELEM", "SPECIAL")
TEX_MAX = 3000          # longer answers are shown verbatim only
RAW_MAX = 6000          # verbatim text is cut here (the cache has it all)


def status_of(v):
    return v if isinstance(v, str) else (v or {}).get("status", "?")


def newest_results(bench):
    fs = sorted(glob.glob(os.path.join(bench, "results", "*.json")), key=os.path.getmtime)
    if not fs:
        sys.exit("no results/*.json: run goursat_stress.py first")
    return json.load(open(fs[-1]))


def answers(bench, methods):
    """(method, id) -> cache record (the last one written: the newest binary)"""
    out = {}
    path = os.path.join(bench, "results", "cache.jsonl")
    for ln in open(path):
        try:
            r = json.loads(ln)
        except ValueError:
            continue
        k = r.get("ckey", "").split("|")
        if len(k) >= 3 and k[0] == "mathilda" and k[1] in methods:
            out[(k[1], k[2])] = r
    return out


def verbatim(s):
    cut = s if len(s) <= RAW_MAX else s[:RAW_MAX] + f" ... [{len(s) - RAW_MAX} more characters]"
    # hard-wrapped: listings cannot break inside one long token (a 10^5-digit integer)
    cut = "\n".join(cut[i:i + 140] for i in range(0, len(cut), 140))
    return ("\\begin{lstlisting}\n" + cut.replace("\\end{lstlisting}", "\\end {lstlisting}") +
            "\n\\end{lstlisting}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", required=True)
    ap.add_argument("--methods", default="GeneralisedGoursat,Trager")
    ap.add_argument("--out", default=None)
    ap.add_argument("--no-pdf", action="store_true")
    a = ap.parse_args()
    methods = a.methods.split(",")
    out = a.out or os.path.join(a.bench, "report")
    os.makedirs(out, exist_ok=True)

    cases = GS.load_cases(os.path.join(a.bench, "corpus.tsv"))
    res = newest_results(a.bench)
    runs = {r["id"]: r["runs"] for r in res["rows"]}
    ans = answers(a.bench, methods)
    lab = {m: f"Mathilda:{m}" for m in methods}

    # TeX for every integrand and every short enough answer
    items = {}
    for c in cases:
        items[f"f|{c['id']}"] = c["f"]
        for m in methods:
            r = ans.get((m, c["id"]))
            if r and r.get("text") and len(r["text"]) <= TEX_MAX:
                items[f"{m}|{c['id']}"] = r["text"]
    print(f"typesetting {len(items)} expressions ...", flush=True)
    tex = tex_batch(items, jobs=8)

    fams = list(dict.fromkeys(c["family"] for c in cases))
    by_fam = collections.defaultdict(list)
    for c in cases:
        by_fam[c["family"]].append(c)

    L = [r"\documentclass[9pt]{extarticle}",
         r"\usepackage[a4paper,margin=14mm]{geometry}",
         r"\usepackage{amsmath,amssymb,booktabs,longtable,listings,xcolor,graphicx}",
         r"\usepackage[hidelinks]{hyperref}",
         r"\lstset{basicstyle=\ttfamily\scriptsize,breaklines=true,breakatwhitespace=false,"
         r"columns=fullflexible,keepspaces=true,showstringspaces=false,xleftmargin=1em}",
         r"\setlength{\parindent}{0pt}\setlength{\parskip}{3pt}\setlength{\emergencystretch}{3em}",
         r"\allowdisplaybreaks",
         r"\title{Benchmark 265-goursat-stress: " + " and ".join(tex_escape(m) for m in methods) +
         r", every integral and answer}",
         r"\date{Run " + tex_escape(str(res.get("date", ""))) + "}",
         r"\begin{document}", r"\maketitle",
         "Every integral of the corpus with the answer each method returned. Status as judged "
         "by \\texttt{PossibleZeroQ[D[F, x] - f]} in Mathilda: ELEM (verified elementary), "
         "SPECIAL (verified, special functions), CERT (certified non-elementary), UNEVAL "
         "(returned unevaluated), TIMEOUT (" + str(res.get("timeout", "?")) + " s), UNVERIFIED, "
         "WRONG. Answers longer than " + str(TEX_MAX) + " characters are shown verbatim only.",
         r"\tableofcontents", r"\clearpage", r"\section{Summary}"]

    # summary: per family status counts and agreement
    sts = ["ELEM", "SPECIAL", "CERT", "UNEVAL", "TIMEOUT", "UNVERIFIED", "WRONG", "ERROR"]
    L.append(r"{\scriptsize\begin{longtable}{ll r " + "r" * len(sts) + "}")
    L.append(r"\toprule family & method & n & " + " & ".join(sts) + r" \\ \midrule\endhead")
    tot = {m: collections.Counter() for m in methods}
    for f in fams + ["total"]:
        cs = cases if f == "total" else by_fam[f]
        for i, m in enumerate(methods):
            cnt = collections.Counter(status_of(runs.get(c["id"], {}).get(lab[m])) for c in cs)
            L.append((tex_escape(f) if i == 0 else "") + " & " + tex_escape(m) + f" & {len(cs)} & " +
                     " & ".join(str(cnt.get(s, 0)) for s in sts) + r" \\")
        L.append(r"\midrule" if f != "total" else r"\bottomrule")
    L.append(r"\end{longtable}}")
    if len(methods) == 2:
        m1, m2 = methods
        L.append(r"\subsection*{Agreement (closed = ELEM or SPECIAL)}")
        L.append(r"{\scriptsize\begin{longtable}{l r r r r r}")
        L.append(r"\toprule family & n & both close & only " + tex_escape(m1) + " & only " +
                 tex_escape(m2) + r" & neither \\ \midrule\endhead")
        for f in fams + ["total"]:
            cs = cases if f == "total" else by_fam[f]
            k = collections.Counter()
            for c in cs:
                r = runs.get(c["id"], {})
                a1, a2 = status_of(r.get(lab[m1])) in CLOSED, status_of(r.get(lab[m2])) in CLOSED
                k[(a1, a2)] += 1
            L.append(f"{tex_escape(f)} & {len(cs)} & {k[(True, True)]} & {k[(True, False)]} & "
                     f"{k[(False, True)]} & {k[(False, False)]} " + r"\\")
        L.append(r"\bottomrule\end{longtable}}")

    # every case
    for f in fams:
        L.append(r"\clearpage\section{" + tex_escape(f) + "}")
        for c in by_fam[f]:
            L.append(r"\subsection*{" + tex_escape(c["id"]) + "}")
            L.append(r"\addcontentsline{toc}{subsection}{" + tex_escape(c["id"]) + "}")
            L.append(r"{\small " + tex_escape(c.get("label", "")) + "; expected: " +
                     tex_escape(c["expect"]) + "}")
            t = tex.get(f"f|{c['id']}")
            L.append(integral(t) if t else verbatim(c["f"]))
            for m in methods:
                r = runs.get(c["id"], {}).get(lab[m])
                st = status_of(r)
                rec = ans.get((m, c["id"]), {})
                tm = rec.get("time")
                lf = rec.get("leaf")
                head = (r"\par\medskip\textbf{" + tex_escape(m) + "}: " + tex_escape(st) +
                        (f", {tm:.3g} s" if isinstance(tm, (int, float)) else "") +
                        (f", LeafCount {lf}" if lf else ""))
                L.append(head)
                text = rec.get("text") or ""
                if st in ("TIMEOUT", "ERROR", "UNEVAL") or not text:
                    if rec.get("detail"):
                        L.append(r"\par{\scriptsize " + tex_escape(rec["detail"][:300]) + "}")
                    continue
                tt = tex.get(f"{m}|{c['id']}")
                if tt:
                    L.append(display("F", tt))
                else:
                    L.append(verbatim(text))
    L.append(r"\end{document}")
    texfile = os.path.join(out, "goursat_vs_trager.tex")
    open(texfile, "w").write("\n".join(L) + "\n")
    print("wrote", texfile, flush=True)
    if not a.no_pdf:
        subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error",
                        os.path.basename(texfile)], cwd=out, capture_output=True)
        pdf = texfile[:-4] + ".pdf"
        print("pdf:", pdf if os.path.exists(pdf) else "FAILED (see the .log)")


if __name__ == "__main__":
    main()
