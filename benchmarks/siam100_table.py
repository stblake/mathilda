#!/usr/bin/env python3
"""
siam100_table.py -- score the SIAM 100-Digit Challenge experiments (166-195)
and their holdouts (196-) against the TRUE answers, in every system.

    # run every system on 166-175, save raw output, write a table
    python3 benchmarks/siam100_table.py --run --only 166-175 \\
        --save /tmp/s100 --md /tmp/s100/SIAM100.md
    # only some systems (the others are read from --from if given)
    python3 benchmarks/siam100_table.py --run --systems mathilda --only 166 \\
        --from /tmp/s100 --save /tmp/s100
    # re-classify saved output
    python3 benchmarks/siam100_table.py --from /tmp/s100 -v

TRUTH.  Each experiment folder holds truth.tsv, one case per line:

    <label> TAB <true value, >= 55 significant digits> TAB <need> TAB <source>

`need` is the number of correct significant digits that makes the case
"correct" (the goal for that tier: typically 10 at machine precision, 45 at
50 digits; lower where the problem itself is ill-conditioned).  `source` is the
provenance (two independent methods, or "consensus: ...").  Lines starting
with # are comments.

SCORING.  Each system's VALUE line {m, e} is read as an exact decimal; its
correct digits are

    digits = floor(-log10(|v - t| / |t|)),  capped at the digits v carries
                                             and at the tier cap (15 / 50)

which is the book's notion of correct digits without its truncation
artefact (0.32336743159999 is ten correct digits of 0.3233674316, not seven).
A case is

    correct   digits >= need
    inexact   an answer with 1 <= digits < need       (counts as not correct)
    wrong     an answer with 0 correct digits, or a complex/malformed one
    none      no answer: timeout, unevaluated, failed, SKIP, or no line

THE TARGET, per case, for Mathilda:

    DIGITS    Mathilda's digits >= the best other system's digits
    SPEED     where Mathilda and at least one other system are correct,
              Mathilda's time < 1.5x the fastest other correct system
              (run_all.py's SLOWER threshold)

A case meeting both is OK; otherwise it is FEWER-DIGITS and/or SLOWER, and
its row says by how much.
"""

import argparse
import json
import math
import os
import re
import sys
from decimal import Decimal, getcontext

getcontext().prec = 120
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_all  # noqa: E402  (discover, run_one, SYSTEMS, PRETTY)

SYSTEMS = ("mathilda", "wolfram", "python", "maxima", "fricas")
PRETTY = {"mathilda": "Mathilda", "wolfram": "WL", "python": "Py",
          "maxima": "Maxima", "fricas": "FriCAS"}
SLOWER_AT = run_all.SLOWER_AT
SIG_RE = re.compile(r"^\{\s*(-?\d+)\s*,\s*(-?\d+)\s*\}$")
EXP_FIRST, EXP_LAST = 166, 999


def tier_cap(label):
    """15 digits for a machine row, the working-precision digits otherwise."""
    m = re.search(r"\bwp(\d+)\b", label)
    return int(m.group(1)) if m else 15


def read_truth(d):
    out = {}
    path = os.path.join(d, "truth.tsv")
    if not os.path.exists(path):
        return out
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln.strip() or ln.lstrip().startswith("#"):
                continue
            parts = ln.split("\t")
            if len(parts) < 3:
                raise SystemExit("%s: malformed truth line: %r" % (path, ln))
            label, val, need = parts[0].strip(), parts[1].strip(), parts[2].strip()
            src = parts[3].strip() if len(parts) > 3 else ""
            out[label] = {"t": Decimal(val), "need": int(need), "src": src}
    return out


def sig_value(s):
    """'{m, e}' -> (Decimal value, digits carried) or None."""
    m = SIG_RE.match(s.strip())
    if not m:
        return None
    mant, e = int(m.group(1)), int(m.group(2))
    if mant == 0:
        return Decimal(0), 1
    nd = len(str(abs(mant)))
    return Decimal(mant).scaleb(e - nd + 1), nd


def digits_correct(v, nd, t, cap):
    if t == 0:
        return min(cap, nd) if v == 0 else 0
    if v == t:
        return min(cap, nd)
    rel = abs(v - t) / abs(t)
    dg = int(math.floor(-float(rel.log10())))
    return max(0, min(dg, nd, cap))


def score(run, label, truth):
    """-> dict(cls, digits, ms, why) for one system on one case."""
    if run is None:
        return {"cls": "none", "digits": None, "ms": None, "why": "not run"}
    if label in run.get("skip", {}):
        return {"cls": "none", "digits": None, "ms": None,
                "why": "SKIP " + run["skip"][label]}
    val = run.get("value", {}).get(label)
    chk = run.get("check", {}).get(label)
    if val is None:
        why = (chk or "no output").strip('"')
        if run.get("timed_out") and chk is None:
            why = "file timeout"
        return {"cls": "none", "digits": None, "ms": None, "why": why}
    sv = sig_value(val)
    if sv is None:
        return {"cls": "wrong", "digits": 0, "ms": None, "why": "malformed"}
    v, nd = sv
    dg = digits_correct(v, nd, truth["t"], tier_cap(label))
    ms = run.get("bench", {}).get(label)
    if dg >= truth["need"]:
        cls = "correct"
    elif dg >= 1:
        cls = "inexact"
    else:
        cls = "wrong"
    return {"cls": cls, "digits": dg, "ms": ms, "why": ""}


def verdict(row):
    """Mathilda's target status for one case."""
    m = row["mathilda"]
    others = [row[s] for s in SYSTEMS if s != "mathilda"]
    best_other = max([o["digits"] or 0 for o in others] + [0])
    md = m["digits"] or 0
    flags = []
    if md < best_other:
        flags.append("FEWER-DIGITS(%d<%d)" % (md, best_other))
    if m["cls"] == "correct":
        times = [o["ms"] for o in others if o["cls"] == "correct" and o["ms"]]
        if times and m["ms"] is not None:
            fastest = min(times)
            if m["ms"] >= SLOWER_AT * fastest:
                flags.append("SLOWER(%.1fx)" % (m["ms"] / fastest))
    if not flags:
        if m["cls"] != "correct" and best_other == 0:
            return "OK(nobody)"
        return "OK"
    return " ".join(flags)


def parse_only(spec):
    if not spec:
        return None
    nums = set()
    for part in spec.split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-")
            nums.update(range(int(a), int(b) + 1))
        elif part:
            nums.add(int(part))
    return nums


def experiments(only):
    out = []
    for e in run_all.discover():
        n = int(e["num"])
        if n < EXP_FIRST or n > EXP_LAST or "siam100" not in e["slug"]:
            continue
        if only and n not in only:
            continue
        out.append(e)
    return out


def fmt_cell(sc):
    if sc["cls"] == "none":
        return "—"
    ms = sc["ms"]
    t = ("%.0f" % ms if ms >= 10 else "%.2g" % ms) if ms is not None else "?"
    mark = {"correct": "", "inexact": "~", "wrong": "✗"}[sc["cls"]]
    return "%s%d d / %s ms" % (mark, sc["digits"], t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--systems", default=",".join(SYSTEMS))
    ap.add_argument("--only", default="")
    ap.add_argument("--from", dest="src", default="")
    ap.add_argument("--save", default="")
    ap.add_argument("--md", default="")
    ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument("-v", action="store_true")
    a = ap.parse_args()

    only = parse_only(a.only)
    exps = experiments(only)
    systems = [s for s in a.systems.split(",") if s]
    runs = {}                         # name -> system -> run dict
    if a.src:
        for e in exps:
            p = os.path.join(a.src, e["name"] + ".json")
            if os.path.exists(p):
                with open(p) as fh:
                    runs[e["name"]] = json.load(fh)
    if a.run:
        for e in exps:
            r = runs.setdefault(e["name"], {})
            for s in systems:
                to = max(a.timeout, e.get("timeout") or 0)
                print("  %-40s %-8s ..." % (e["name"], s), file=sys.stderr,
                      flush=True)
                r[s] = run_all.run_one(s, e, to)
            if a.save:
                os.makedirs(a.save, exist_ok=True)
                with open(os.path.join(a.save, e["name"] + ".json"), "w") as fh:
                    json.dump(r, fh, indent=1)

    lines, tally = [], {}
    head = "| case | need | " + " | ".join(PRETTY[s] for s in SYSTEMS) + " | verdict |"
    sep = "|" + "---|" * (len(SYSTEMS) + 3)
    summary = []
    for e in exps:
        truth = read_truth(e["dir"])
        r = runs.get(e["name"], {})
        labels = list(truth)
        extra = set()
        for s in SYSTEMS:
            rr = r.get(s) or {}
            extra |= set(rr.get("check", {})) | set(rr.get("skip", {}))
        missing = sorted(extra - set(truth))
        if missing:
            print("WARNING %s: labels with no truth: %s" % (e["name"], missing),
                  file=sys.stderr)
        lines.append("\n### %s\n\n%s\n%s" % (e["name"], head, sep))
        ok = 0
        for lb in labels:
            row = {s: score(r.get(s), lb, truth[lb]) for s in SYSTEMS}
            vd = verdict(row)
            ok += vd.startswith("OK")
            for s in SYSTEMS:
                tally.setdefault(s, {"correct": 0, "inexact": 0, "wrong": 0,
                                     "none": 0})[row[s]["cls"]] += 1
            lines.append("| %s | %d | %s | %s |" % (
                lb, truth[lb]["need"],
                " | ".join(fmt_cell(row[s]) for s in SYSTEMS), vd))
            if a.v and not vd.startswith("OK"):
                print("%-28s %s   [%s]" % (lb, vd, "; ".join(
                    "%s:%s" % (s, row[s]["why"] or row[s]["cls"])
                    for s in SYSTEMS)), file=sys.stderr)
        summary.append("| %s | %d | %d |" % (e["name"], len(labels), ok))

    out = ["# SIAM 100-digit scoring", "",
           "Cell: correct digits / ms (`~` below need, `✗` wrong, — none).", "",
           "| experiment | cases | Mathilda on target |", "|---|---|---|"]
    out += summary
    out += ["", "| system | correct | inexact | wrong | none |", "|---|---|---|---|---|"]
    for s in SYSTEMS:
        t = tally.get(s, {"correct": 0, "inexact": 0, "wrong": 0, "none": 0})
        out.append("| %s | %d | %d | %d | %d |" % (PRETTY[s], t["correct"],
                                                   t["inexact"], t["wrong"], t["none"]))
    out += lines
    text = "\n".join(out) + "\n"
    if a.md:
        os.makedirs(os.path.dirname(os.path.abspath(a.md)), exist_ok=True)
        with open(a.md, "w") as fh:
            fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
