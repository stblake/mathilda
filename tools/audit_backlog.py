#!/usr/bin/env python3
"""
audit_backlog.py -- the live backlog of the 2026-11-14 builtin audit.

The audit (review/report.html, one write-up per builtin in review/issues/)
confirmed 2,475 defects against the binary of that day.  The tree has moved
since, so before anything is fixed every finding is re-run against the
CURRENT ./Mathilda and classified:

  CRASH            the process died on a signal
  HANG             killed by the time limit or the RSS cap
  SAME             output identical to the audit's "Actual" (still broken)
  CHANGED          output differs from "Actual"; needs adjudication against
                   "Expected" (may be fixed, may be differently broken)
  LIKELY-FIXED     CHANGED, and every result line occurs in "Expected"

Each finding runs in its own process, one statement per line, printed as
  @@<i> <InputForm of the result>
with stderr merged so messages interleave exactly as the audit's mrun
output did.  Every child is watched for wall time AND resident memory:
TimeConstrained cannot stop a single GMP call, and an unwatched probe
driver once reached 47 GB (see check_machine_boundary.py).

Usage:
  tools/audit_backlog.py parse                 # -> review/backlog.json
  tools/audit_backlog.py run [-j4] [--sev=S] [--builtin=B] [--ids=a#1,b#2]
                             [--resume]       # skip ids already in the status
                                               # -> review/backlog_status.json/.md
  tools/audit_backlog.py report                # regenerate the .md only

Environment: AB_TIMEOUT (s, default 30), AB_MEM_CAP_MB (default 2048).
"""

import concurrent.futures as cf
import json
import os
import re
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ISSUES = os.path.join(ROOT, "review", "issues")
BACKLOG = os.path.join(ROOT, "review", "backlog.json")
STATUS_JSON = os.path.join(ROOT, "review", "backlog_status.json")
STATUS_MD = os.path.join(ROOT, "review", "backlog_status.md")
BIN = os.path.join(ROOT, "Mathilda")
TIMEOUT = float(os.environ.get("AB_TIMEOUT", "30"))
MEM_CAP_MB = float(os.environ.get("AB_MEM_CAP_MB", "2048"))
POLL_SECONDS = 0.05
SEV_RANK = {"critical": 0, "high": 1, "medium": 2, "low": 3}
AUDITOR_PREFIX = re.compile(r"/Users/[^/]+/Downloads/[^/]+/[^/]+/")

# Tier A: the core architecture (CLAUDE.md: parser, matcher, symbol table,
# evaluator, the simple arithmetic heads) and the builtins the audit itself
# ranked first (audit-summary-batches-1-2 / 3-5 name them).
TIER_A = set("""
Power Times Integrate Log Plus Table D Solve Sort Cot Position Refine Root Module
NIntegrate Part Reduce Expand Piecewise Apart Cancel Flatten With Replace NSum
Set SetDelayed Equal SameQ UnsameQ Map Apply Function Slot SlotSequence Block
ReplaceAll ReplaceRepeated Rule RuleDelayed List Length Head Simplify Factor
Together Series Limit Sum Product Divide Subtract Minus Sqrt Exp Sin Cos Tan
If Which Do For While Return Break Continue Catch Throw Scan MatchQ Cases
Select Join Range Rational Complex N Abs Less Greater LessEqual GreaterEqual
Unequal And Or Not Clear ClearAll Protect Unprotect Attributes SetAttributes
Blank BlankSequence BlankNullSequence Pattern Condition PatternSequence
Optional Hold HoldForm Unevaluated Evaluate Nest NestList Fold FoldList
FixedPoint Thread MapThread Total Mean Max Min Floor Ceiling Round Mod
""".split())
TIER_B = set("""
ReplacePart Factor Zeta LinearSolve ComplexExpand If Join Gamma Derivative
Eigenvalues Inverse Norm PolyLog Floor Ceiling Det NSolve MatrixRank PowerMod
Erf BesselJ BesselK EllipticE EllipticF EllipticK LegendreP Residue Partition
LeastSquares Information FindInstance Delete Insert Fit Quantile LegendreQ
NMinimize NMaximize FindMinimum Eigensystem Plot ParametricPlot Coefficient
CoefficientList Collect Exponent PolynomialQuotient PolynomialRemainder
PolynomialGCD Dot Transpose Take Drop First Last Rest Most Append Prepend
Reverse Union Intersection Complement Count Tally DeleteDuplicates StringJoin
StringReplace StringSplit ToString ToExpression Print Timing Variables
FullSimplify TrigExpand TrigReduce Numerator Denominator Re Im Conjugate Arg
""".split())


# --------------------------------------------------------------------- parse

FINDING = re.compile(
    r"^## (\d+)\. (.*?)\n\n\*\*Kind:\*\* (\w+) · \*\*Severity:\*\* (\w+)"
    r"(?: · \*\*Likely cause:\*\* `(.*?)`)?(.*?)(?=^## \d+\. |\Z)",
    re.S | re.M)


def _block(body, label):
    m = re.search(r"\*\*" + label + r":?\*\*[^\n]*\n```[a-z]*\n(.*?)```",
                  body, re.S)
    return m.group(1) if m else None


def _para(body, label):
    m = re.search(r"\*\*" + label + r":\*\*\s*(.*?)(?=\n\n\*\*|\Z)", body, re.S)
    return m.group(1).strip() if m else None


def split_statements(src):
    """Group lines into complete statements: brackets balanced outside
    strings and comments.  Pure comment lines are dropped."""
    out, cur, depth, in_str, in_cmt = [], [], 0, False, 0
    for line in src.split("\n"):
        if not cur and not line.strip():
            continue
        cur.append(line)
        i = 0
        while i < len(line):
            c = line[i]
            two = line[i:i + 2]
            if in_cmt:
                if two == "*)":
                    in_cmt -= 1; i += 1
                elif two == "(*":
                    in_cmt += 1; i += 1
            elif in_str:
                if c == "\\":
                    i += 1
                elif c == '"':
                    in_str = False
            elif two == "(*":
                in_cmt += 1; i += 1
            elif c == '"':
                in_str = True
            elif c in "([{":
                depth += 1
            elif c in ")]}":
                depth -= 1
            i += 1
        if depth <= 0 and not in_str and not in_cmt:
            stmt = "\n".join(cur).strip()
            cur, depth = [], 0
            if re.fullmatch(r"\(\*.*\*\)", stmt, re.S):
                continue
            if stmt:
                out.append(stmt)
    if cur:
        out.append("\n".join(cur).strip())
    return out


def parse():
    rows = []
    for fn in sorted(os.listdir(ISSUES)):
        if not fn.endswith(".md"):
            continue
        b = fn[:-3]
        text = open(os.path.join(ISSUES, fn), encoding="utf-8").read()
        cat = re.search(r"\*\*Category:\*\* (\S+)", text)
        reg = re.search(r"\*\*Registered in:\*\* `(.*?)`", text)
        for m in FINDING.finditer(text):
            n, title, kind, sev, cause, body = m.groups()
            rep = _block(body, "Reproduce")
            rows.append({
                "id": f"{b}#{n}", "builtin": b, "n": int(n), "title": title,
                "kind": kind, "severity": sev,
                "category": cat.group(1) if cat else None,
                "registered": reg.group(1) if reg else None,
                "cause": AUDITOR_PREFIX.sub("", cause) if cause else None,
                "statements": split_statements(rep) if rep else [],
                "actual": _block(body, "Actual"),
                "expected": _para(body, "Expected"),
                "tier": "A" if b in TIER_A else "B" if b in TIER_B else "C",
            })
    rows.sort(key=lambda r: (SEV_RANK[r["severity"]], r["tier"], r["id"]))
    json.dump(rows, open(BACKLOG, "w"), indent=1)
    return rows


# ----------------------------------------------------------------------- run

def rss_mb(pid):
    try:
        out = subprocess.run(["ps", "-o", "rss=", "-p", str(pid)],
                             capture_output=True, text=True).stdout
        return int(out.strip() or 0) / 1024.0
    except (ValueError, OSError):
        return 0.0


def watch(p):
    """Wait for p, killing it past TIMEOUT or MEM_CAP_MB."""
    t0 = time.monotonic()
    k = 0
    while p.poll() is None:
        k += 1
        if k % 4 == 0 and rss_mb(p.pid) > MEM_CAP_MB:
            p.kill(); p.wait()
            return f"MEMCAP>{MEM_CAP_MB:.0f}MB"
        if time.monotonic() - t0 > TIMEOUT:
            p.kill(); p.wait()
            return "TIMEOUT"
        time.sleep(POLL_SECONDS)
    return None


def script_for(stmts):
    lines = []
    for i, s in enumerate(stmts, 1):
        body = s.rstrip()
        if body.endswith(";"):
            # A suppressed statement: run it, report Null like the REPL did.
            lines.append(f"{body[:-1]};")
            lines.append(f'Print["@@{i} Null"];')
        else:
            lines.append(f'Print["@@{i} ", ToString[InputForm[{body}]]];')
    return "\n".join(lines) + "\n"


def run_one(row):
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as fh:
        fh.write(script_for(row["statements"]))
        path = fh.name
    outf = tempfile.TemporaryFile("w+b")
    t0 = time.monotonic()
    try:
        p = subprocess.Popen([BIN, "-file", path], stdout=outf,
                             stderr=subprocess.STDOUT,
                             cwd=tempfile.gettempdir())
        err = watch(p)
        outf.seek(0)
        out = outf.read().decode("utf-8", errors="replace")
        rc = p.returncode
    finally:
        os.unlink(path)
        outf.close()
    return out, rc, err, time.monotonic() - t0


def _norm(s):
    return re.sub(r"\s+", " ", s).strip()


def _tight(s):
    return re.sub(r"[\s*]", "", s)


def actual_lines(actual):
    """The audit's output, reduced to the comparable lines."""
    keep = []
    for line in (actual or "").split("\n"):
        if re.match(r"^\[\d+\] In:", line) or not line.strip():
            continue
        m = re.match(r"^\[\d+\] Out= ?(.*)$", line)
        keep.append(("OUT", _norm(m.group(1))) if m else ("MSG", _norm(line)))
    return keep


def current_lines(out):
    keep = []
    for line in out.split("\n"):
        if not line.strip():
            continue
        m = re.match(r"^@@\d+ ?(.*)$", line)
        keep.append(("OUT", _norm(m.group(1))) if m else ("MSG", _norm(line)))
    return keep


def classify(row, out, rc, err):
    if err:
        return "HANG", err
    if rc is not None and rc < 0:
        return "CRASH", f"signal {-rc}"
    cur = current_lines(out)
    act = actual_lines(row["actual"])
    if "<<CRASH" in (row["actual"] or "") or "<<TIMEOUT" in (row["actual"] or ""):
        return "CHANGED", "was crash/timeout, now terminates"
    if cur == act or [_tight(v) for _, v in cur] == [_tight(v) for _, v in act]:
        return "SAME", ""
    exp = _tight(row["expected"] or "")
    outs = [_tight(v) for k, v in cur if k == "OUT"]
    if outs and exp and all(o and o in exp for o in outs):
        return "LIKELY-FIXED", ""
    return "CHANGED", ""


def run(args):
    rows = json.load(open(BACKLOG)) if os.path.exists(BACKLOG) else parse()
    jobs = 4
    sev = builtin = ids = None
    resume = False
    for a in args:
        if a.startswith("-j"):
            jobs = min(4, int(a[2:]))
        elif a.startswith("--sev="):
            sev = set(a[6:].split(","))
        elif a.startswith("--builtin="):
            builtin = set(a[10:].split(","))
        elif a.startswith("--ids="):
            ids = set(a[6:].split(","))
        elif a == "--resume":
            resume = True
    sel = [r for r in rows
           if (not sev or r["severity"] in sev)
           and (not builtin or r["builtin"] in builtin)
           and (not ids or r["id"] in ids)]
    status = json.load(open(STATUS_JSON)) if os.path.exists(STATUS_JSON) else {}
    if resume:
        sel = [r for r in sel if r["id"] not in status]
    done = 0
    with cf.ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = {ex.submit(run_one, r): r for r in sel}
        for f in cf.as_completed(futs):
            r = futs[f]
            try:
                out, rc, err, secs = f.result()
                st, note = classify(r, out, rc, err)
            except Exception as ex:          # a harness fault, not a verdict
                out, rc, secs, st, note = "", None, 0.0, "HARNESS-ERROR", repr(ex)
            status[r["id"]] = {"status": st, "note": note, "rc": rc,
                               "secs": round(secs, 2), "output": out[-4000:]}
            done += 1
            if done % 50 == 0:           # checkpoint: a later failure keeps these
                json.dump(status, open(STATUS_JSON, "w"), indent=1)
            if done % 50 == 0 or st in ("CRASH", "HANG"):
                print(f"[{done}/{len(sel)}] {r['id']}: {st} {note}",
                      flush=True)
    json.dump(status, open(STATUS_JSON, "w"), indent=1)
    report(rows, status)


def report(rows=None, status=None):
    rows = rows or json.load(open(BACKLOG))
    status = status or json.load(open(STATUS_JSON))
    from collections import Counter
    c = Counter()
    for r in rows:
        s = status.get(r["id"], {}).get("status", "NOT-RUN")
        c[(r["severity"], s)] += 1
    states = ["CRASH", "HANG", "SAME", "CHANGED", "LIKELY-FIXED", "NOT-RUN"]
    md = ["# Audit backlog status", "",
          f"*Generated {time.strftime('%Y-%m-%d %H:%M')} by "
          "`tools/audit_backlog.py` against the current `./Mathilda`.*", "",
          "| severity | " + " | ".join(states) + " |",
          "|---|" + "---|" * len(states)]
    for sv in SEV_RANK:
        md.append(f"| {sv} | " + " | ".join(str(c[(sv, s)]) for s in states)
                  + " |")
    md += ["", "## Open findings (crash, hang, same), by priority", "",
           "| id | sev | tier | kind | status | title |", "|---|---|---|---|---|---|"]
    for r in rows:
        s = status.get(r["id"], {}).get("status", "NOT-RUN")
        if s in ("CRASH", "HANG", "SAME"):
            md.append(f"| {r['id']} | {r['severity']} | {r['tier']} | "
                      f"{r['kind']} | {s} | {r['title'][:90].replace('|', '/')} |")
    open(STATUS_MD, "w").write("\n".join(md) + "\n")
    print("\n".join(md[:12]))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "parse":
        rs = parse()
        print(f"{len(rs)} findings -> {BACKLOG}")
    elif cmd == "run":
        run(sys.argv[2:])
    elif cmd == "report":
        report()
    else:
        sys.exit(__doc__)
