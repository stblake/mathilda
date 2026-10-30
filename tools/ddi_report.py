#!/usr/bin/env python3
"""
ddi_report.py -- render tests/data/ddi_vol1_report.json (written by
tools/ddi_stress.py) as the data sections of ALGEBRAIC_INTEGRAL_PERF.md.

Every number in the generated tables is computed from the JSON; the prose
(analysis, gap list) lives in a hand-written header/footer which this tool
splices around the generated block:

    python3 tools/ddi_report.py [--json J] [--md ALGEBRAIC_INTEGRAL_PERF.md]

The generated block sits between the markers
    <!-- BEGIN GENERATED -->  and  <!-- END GENERATED -->
and is replaced in place (the file is created with an empty frame if absent).
"""
import argparse
import json
import os
import re
import statistics

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BEGIN, END = "<!-- BEGIN GENERATED -->", "<!-- END GENERATED -->"
CLASSES = ("AGREE-ELEM", "MATHILDA-ONLY", "AGREE-NONELEM", "NONELEM-UNCONFIRMED",
           "FRICAS-ONLY", "NEITHER", "TIMEOUT", "ERROR", "WRONG", "BRANCH", "PZQ-FN", "UNVERIFIED",
           "DISAGREE", "BOOK-CONFLICT")
SHORT = {"AGREE-ELEM": "agree elem", "MATHILDA-ONLY": "Mathilda only",
         "AGREE-NONELEM": "agree non-elem", "NONELEM-UNCONFIRMED": "non-elem unconf.",
         "FRICAS-ONLY": "FriCAS only", "NEITHER": "neither", "TIMEOUT": "timeout",
         "ERROR": "error", "WRONG": "WRONG", "BRANCH": "wrong branch on part of R",
         "PZQ-FN": "PossibleZeroQ false negative", "UNVERIFIED": "unverified", "DISAGREE": "DISAGREE",
         "BOOK-CONFLICT": "book conflict"}


def cell(s, n=110):
    s = "" if s is None else str(s)
    s = s if len(s) <= n else s[:n] + "…"
    return "`" + s.replace("`", "'").replace("|", "\\|").replace("\n", " ") + "`" if s else ""


def esc(s, n=110):
    """Plain text with every markdown-active character escaped, so a long
    integrand WRAPS inside a table cell (a code span cannot break in LaTeX)."""
    s = "" if s is None else str(s)
    s = s if len(s) <= n else s[:n] + "…"
    s = s.replace("\n", " ")
    return re.sub(r"([\\`*_{}\[\]()#+!|^~<>$&-])", r"\\\1", s)


def plain(s, n=110):
    s = "" if s is None else str(s)
    s = s if len(s) <= n else s[:n] + "…"
    return s.replace("|", "\\|").replace("\n", " ")


def tstats(xs):
    xs = [v for v in xs if v is not None]
    if not xs:
        return ["-"] * 4
    return [str(len(xs)), f"{statistics.median(xs):.3f}", f"{statistics.mean(xs):.3f}",
            f"{max(xs):.2f}"]


def solved(m):
    return m["status"] == "ELEM"


def decided(m):
    return m["status"] in ("ELEM", "NONELEM")


def generate(d):
    cases, methods = d["cases"], d["methods"]
    L = []
    A = L.append
    fam_order = []
    for c in cases:
        if c["family"] not in fam_order:
            fam_order.append(c["family"])
    nsym = sum(1 for c in cases if c["params"])
    A(f"Run of {d['generated']}: **{len(cases)} integrands** "
      f"({nsym} with symbolic parameters, {len(cases) - nsym} over Q), "
      f"per-method timeout {d['timeout']} s, FriCAS timeout {d['fricas_timeout']} s, "
      f"wall time {d['elapsed'] / 60:.0f} min" +
      ("".join(f"; {', '.join(g['methods'])} re-run on {g['cases']} cases "
               f"({g['elapsed'] / 60:.0f} min)" for g in d.get("merged", []))) + ".\n")

    # --- headline -----------------------------------------------------------
    A("### Headline\n")
    hdr = ["", *methods, "FriCAS"]
    A("| " + " | ".join(hdr) + " |")
    A("|" + "---|" * len(hdr))

    def frow(name, fn_m, fn_f):
        A("| " + " | ".join([name] + [str(sum(1 for c in cases if fn_m(c["methods"][m])))
                                     for m in methods] +
                            [str(sum(1 for c in cases if fn_f(c["fricas"])))]) + " |")
    frow("verified antiderivative", solved, lambda f: f["status"] == "ELEM")
    frow("proved non-elementary", lambda m: m["status"] == "NONELEM",
         lambda f: f["status"] == "NONELEM")
    frow("declined", lambda m: m["status"] == "DECLINE", lambda f: False)
    frow("timeout", lambda m: m["status"] == "TIMEOUT", lambda f: f["status"] == "TIMEOUT")
    frow("error / crash", lambda m: m["status"] == "ERROR", lambda f: f["status"] == "ERROR")
    frow("wrong answer", lambda m: m["status"] == "WRONG", lambda f: f["status"] == "BAD")
    frow("right on part of the real line only (branch)", lambda m: m["status"] == "BRANCH",
         lambda f: False)
    frow("correct, but PossibleZeroQ cannot confirm it (PZQ-FN)",
         lambda m: m["status"] == "PZQ-FN", lambda f: False)
    frow("answer not verifiable in budget", lambda m: m["status"] == "UNVERIFIED",
         lambda f: False)
    A("")
    e_book = [c for c in cases if c["expect"] == "elem"]
    if e_book:
        A(f"Of the **{len(e_book)}** integrands the book states are elementary: " +
          ", ".join(f"{m} solves **{sum(1 for c in e_book if solved(c['methods'][m]))}**"
                    for m in methods) +
          f", FriCAS **{sum(1 for c in e_book if c['fricas']['status'] == 'ELEM')}**, "
          f"either Mathilda method **"
          f"{sum(1 for c in e_book if any(solved(c['methods'][m]) for m in methods))}**.\n")
    for sub, pred in (("over Q (no parameters)", lambda c: not c["params"]),
                      ("with symbolic parameters", lambda c: c["params"])):
        cs = [c for c in cases if pred(c)]
        if cs:
            A(f"- {sub} ({len(cs)}): " + "; ".join(
                f"{m} solved {sum(1 for c in cs if solved(c['methods'][m]))}, "
                f"decided {sum(1 for c in cs if decided(c['methods'][m]))}"
                for m in methods) +
              f"; FriCAS solved {sum(1 for c in cs if c['fricas']['status'] == 'ELEM')}.")
    A("")

    # --- classification ------------------------------------------------------
    A("### Classification against FriCAS and the book\n")
    A("| class | " + " | ".join(methods) + " |")
    A("|---|" + "---:|" * len(methods))
    for k in CLASSES:
        row = [str(sum(1 for c in cases if c["methods"][m]["class"] == k)) for m in methods]
        if any(r != "0" for r in row):
            A(f"| {SHORT[k]} | " + " | ".join(row) + " |")
    A("")

    # --- per chapter ---------------------------------------------------------
    A("### Per chapter\n")
    A("Cells are *verified antiderivatives / proved non-elementary / declined / "
      "timeout+error+wrong*.\n")
    A("| chapter | n | book elem | " + " | ".join(methods) + " | FriCAS (elem/non) |")
    A("|---|---:|---:|" + "---|" * (len(methods) + 1))
    for fam in fam_order:
        cs = [c for c in cases if c["family"] == fam]
        cells = []
        for m in methods:
            s = [c["methods"][m]["status"] for c in cs]
            cells.append(f"{s.count('ELEM')} / {s.count('NONELEM')} / {s.count('DECLINE')} / "
                         f"{s.count('TIMEOUT') + s.count('ERROR') + s.count('WRONG')}")
        fs = [c["fricas"]["status"] for c in cs]
        A(f"| {fam} | {len(cs)} | {sum(1 for c in cs if c['expect'] == 'elem')} | " +
          " | ".join(cells) + f" | {fs.count('ELEM')} / {fs.count('NONELEM')} |")
    A("")

    # --- timing --------------------------------------------------------------
    A("### Timing (seconds, cases each system answered — antiderivative or proof)\n")
    A("| system | answered | median | mean | max |")
    A("|---|---:|---:|---:|---:|")
    for m in methods:
        A(f"| {m} | " + " | ".join(tstats([c["methods"][m].get("time") for c in cases
                                            if decided(c["methods"][m])])) + " |")
    A("| FriCAS | " + " | ".join(tstats([c["fricas"].get("time") for c in cases
                                         if c["fricas"]["status"] in ("ELEM", "NONELEM")])) + " |")
    A("")
    both = [c for c in cases if all(solved(c["methods"][m]) for m in methods)]
    if both:
        lv = {m: [c["methods"][m]["leaves"] for c in both if c["methods"][m].get("leaves")]
              for m in methods}
        A(f"On the {len(both)} integrands both methods solve, median LeafCount of the "
          "antiderivative: " + ", ".join(f"{m} {statistics.median(lv[m]):.0f}"
                                         for m in methods if lv[m]) + ".\n")

    # --- decline reasons -------------------------------------------------------
    A("### Why each method declines\n")
    for m in methods:
        reasons = {}
        for c in cases:
            r = c["methods"][m]
            if r["status"] == "DECLINE":
                # {"failed", "reason", ...} -> reason; {"reason", data...} -> reason
                strs = re.findall(r'"([^"]*)"', r.get("detail") or "")
                key = (strs[1] if strs and strs[0] == "failed" and len(strs) > 1
                       else strs[0] if strs else (r.get("detail") or "")[:60])
                key = key.split(":")[0][:90]
                reasons[key] = reasons.get(key, 0) + 1
        A(f"**{m}** ({sum(reasons.values())} declines)\n")
        A("| reason | count |")
        A("|---|---:|")
        for k, v in sorted(reasons.items(), key=lambda kv: -kv[1]):
            A(f"| {plain(k)} | {v} |")
        A("")

    # --- listings ---------------------------------------------------------------
    def listing(title, pred, extra=None):
        rows = [(c, m) for c in cases for m in methods if pred(c, c["methods"][m])]
        A(f"### {title} ({len(rows)})\n")
        if not rows:
            A("None.\n")
            return
        A("| id | method | integrand | book ref | detail |")
        A("|-----------|------|------------------|---------|--------------|")
        for c, m in rows:
            r = c["methods"][m]
            det = extra(c, r) if extra else r.get("detail")
            A(f"| {esc(c['id'], 40)} | {'PMT' if m == 'ParallelMixedTower' else m} | "
              f"{esc(c['f'], 90)} | {esc(c['ref'], 40)} | {esc(det, 110)} |")
        A("")

    listing("Wrong answers (failed derivative check)", lambda c, r: r["status"] == "WRONG",
            lambda c, r: f"{r.get('detail') or ''} num={r.get('numeric')}")
    listing("Right on part of the real line only (D[F] = f on one interval, -f on another)",
            lambda c, r: r["status"] == "BRANCH",
            lambda c, r: f"{r.get('detail')}; points {r.get('numeric')}")
    listing("Correct answers PossibleZeroQ cannot confirm (PossibleZeroQ bugs)",
            lambda c, r: r["status"] == "PZQ-FN",
            lambda c, r: f"{r.get('detail')}; F = {r.get('result')}")
    listing("Antiderivatives not verifiable within the budget", lambda c, r: r["status"] == "UNVERIFIED",
            lambda c, r: f"{r.get('detail')}; {r.get('result')}")
    listing("Disagreements (proved non-elementary, but FriCAS or the book has an "
            "antiderivative)", lambda c, r: r["class"] == "DISAGREE",
            lambda c, r: f"{r.get('detail')}; FriCAS {c['fricas']['status']}; book {c['expect'] or '?'}")
    listing("Book conflicts (verified antiderivative where the book says non-elementary)",
            lambda c, r: r["class"] == "BOOK-CONFLICT", lambda c, r: r.get("result"))
    listing("Crashes / errors", lambda c, r: r["status"] == "ERROR")
    listing("Timeouts", lambda c, r: r["status"] == "TIMEOUT",
            lambda c, r: f"FriCAS {c['fricas']['status']}; book {c['expect'] or '?'}")
    listing("Elementary per book or FriCAS, but declined",
            lambda c, r: r["status"] == "DECLINE" and
            (c["expect"] == "elem" or c["fricas"]["status"] == "ELEM"))
    rows = [c for c in cases if c["fricas"]["status"] == "BAD"]
    A(f"### FriCAS answers that fail their own derivative check ({len(rows)})\n")
    for c in rows:
        A(f"- `{c['id']}`: {cell(c['f'], 80)}")
    A("")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=os.path.join(REPO, "tests", "data", "ddi_vol1_report.json"))
    ap.add_argument("--md", default=os.path.join(REPO, "ALGEBRAIC_INTEGRAL_PERF.md"))
    args = ap.parse_args()
    d = json.load(open(args.json))
    gen = generate(d)
    if os.path.exists(args.md):
        text = open(args.md).read()
    else:
        text = f"# Algebraic integration stress test\n\n{BEGIN}\n{END}\n"
    if BEGIN not in text:
        text += f"\n{BEGIN}\n{END}\n"
    head, rest = text.split(BEGIN, 1)
    tail = rest.split(END, 1)[1] if END in rest else ""
    open(args.md, "w").write(head + BEGIN + "\n" + gen + "\n" + END + tail)
    print(f"wrote {args.md}")


if __name__ == "__main__":
    main()
