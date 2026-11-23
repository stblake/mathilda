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
  tools/audit_backlog.py report                # regenerate the .md and the html
  tools/audit_backlog.py pinned [--dry]        # CHANGED rows pinned by a passing
                                               # tests/test_audit_*.c check -> FIXED
  tools/audit_backlog.py html                  # regenerate review/report.html only
                                               # (from review/report.source.html,
                                               # stamped with each finding's status)

Environment: AB_TIMEOUT (s, default 30), AB_MEM_CAP_MB (default 2048),
AB_BIN (binary to check, default ./Mathilda; e.g. a worktree's).
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
# Hand/agent verdicts on CHANGED rows: {"<id>": {"verdict": FIXED |
# STILL-BROKEN | PARTIAL | FORMAT-ONLY | UNCLEAR, "reason": ...}}.  A verdict
# applies only while the row's status is still CHANGED / LIKELY-FIXED.
ADJUDICATION = os.path.join(ROOT, "review", "adjudication.json")
OPEN_STATES = ("CRASH", "HANG", "SAME", "STILL-BROKEN", "PARTIAL",
               "FORMAT-ONLY", "UNCLEAR", "CHANGED", "LIKELY-FIXED")
# AB_BIN: check another tree's binary (a branch before it lands); its
# src/version.h is read alongside it.
BIN = os.path.abspath(os.environ.get("AB_BIN", os.path.join(ROOT, "Mathilda")))
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
    # what this status was measured against (read back by html_report)
    status["_meta"] = {"version": _bin_version(),
                       "time": time.strftime("%Y-%m-%d %H:%M"),
                       "scope": "partial" if (sev or builtin or ids) else "full"}
    json.dump(status, open(STATUS_JSON, "w"), indent=1)
    report(rows, status)


def effective_status(rid, status, adj):
    """The run's classification, refined by an adjudication verdict."""
    s = status.get(rid, {}).get("status", "NOT-RUN")
    if s in ("CHANGED", "LIKELY-FIXED") and rid in adj:
        return adj[rid]["verdict"]
    return s


def report(rows=None, status=None):
    rows = rows or json.load(open(BACKLOG))
    status = status or json.load(open(STATUS_JSON))
    adj = json.load(open(ADJUDICATION)) if os.path.exists(ADJUDICATION) else {}
    from collections import Counter
    c = Counter()
    for r in rows:
        s = effective_status(r["id"], status, adj)
        c[(r["severity"], s)] += 1
    states = ["CRASH", "HANG", "SAME", "STILL-BROKEN", "PARTIAL", "FORMAT-ONLY",
              "UNCLEAR", "CHANGED", "LIKELY-FIXED", "FIXED", "NOT-RUN"]
    md = ["# Audit backlog status", "",
          f"*Generated {time.strftime('%Y-%m-%d %H:%M')} by "
          "`tools/audit_backlog.py` against the current `./Mathilda`.*", "",
          "| severity | " + " | ".join(states) + " |",
          "|---|" + "---|" * len(states)]
    for sv in SEV_RANK:
        md.append(f"| {sv} | " + " | ".join(str(c[(sv, s)]) for s in states)
                  + " |")
    md += ["", "## Open findings, by priority", "",
           "| id | sev | tier | kind | status | title |", "|---|---|---|---|---|---|"]
    for r in rows:
        s = effective_status(r["id"], status, adj)
        if s in OPEN_STATES:
            md.append(f"| {r['id']} | {r['severity']} | {r['tier']} | "
                      f"{r['kind']} | {s} | {r['title'][:90].replace('|', '/')} |")
    open(STATUS_MD, "w").write("\n".join(md) + "\n")
    print("\n".join(md[:12]))
    html_report(rows, status, adj)


def adjudicate_pinned(args):
    """Verdict FIXED for every CHANGED / LIKELY-FIXED row whose id is pinned in
    a tests/test_audit_<category>.c check of the checked tree.

    Each fix is pinned by audit id with checks written from the finding's
    "Expected", so a row that has moved off the audit's "Actual" AND is pinned
    in a passing suite is fixed.  Run this only after the audit suites of that
    tree pass (ctest -R audit_); a row that is not pinned keeps waiting for a
    hand verdict.  --dry lists without writing."""
    dry = "--dry" in args
    tree = os.path.dirname(BIN)
    status = json.load(open(STATUS_JSON))
    adj = json.load(open(ADJUDICATION)) if os.path.exists(ADJUDICATION) else {}
    pinned = {}
    tdir = os.path.join(tree, "tests")
    for fn in sorted(os.listdir(tdir)):
        if fn.startswith("test_audit_") and fn.endswith(".c"):
            for rid in re.findall(r'audit_\w+\(\s*"([^"]+#\d+)"',
                                  open(os.path.join(tdir, fn)).read()):
                pinned.setdefault(rid, fn)
    n = 0
    for rid, st in status.items():
        if rid.startswith("_") or st.get("status") not in ("CHANGED", "LIKELY-FIXED"):
            continue
        if rid in pinned and adj.get(rid, {}).get("verdict") != "FIXED":
            n += 1
            print(f"{rid}: FIXED (pinned in tests/{pinned[rid]})")
            adj[rid] = {"verdict": "FIXED",
                        "reason": f"Output moved off the audit's Actual and the id is "
                                  f"pinned by passing checks in tests/{pinned[rid]} "
                                  f"(v{_bin_version()})."}
    if not dry:
        json.dump(adj, open(ADJUDICATION, "w"), indent=1)
        report()
    print(f"{n} rows adjudicated FIXED from pinned tests"
          + (" (dry run)" if dry else ""))


# --- the annotated audit report ---------------------------------------------
#
# review/report.html is the audit's own write-up with every finding stamped by
# its live status.  It is regenerated from the pristine copy
# review/report.source.html (made once, from the original report.html), so the
# rendering is idempotent: run it as often as the status changes.  The audit
# lists builtins in <details> blocks and each finding in a <div class='f'>, in
# the order of review/issues/<Head>.md -- so the k-th finding of builtin H is
# the backlog row "H#k" (checked: every title matches).

# Display group of each effective status: FIXED is done; CHANGED / LIKELY-FIXED
# / UNCLEAR / FORMAT-ONLY / NOT-RUN still need a human or agent verdict.
HTML_GROUP = {"FIXED": "fixed",
              "CRASH": "open", "HANG": "open", "SAME": "open",
              "STILL-BROKEN": "open", "PARTIAL": "open"}
HTML_SOURCE = os.path.join(ROOT, "review", "report.source.html")
HTML_OUT = os.path.join(ROOT, "review", "report.html")

HTML_STYLE = """
.st{display:inline-block;padding:0 .45em;border-radius:4px;font-size:.75em;font-weight:600;margin-left:.4em;vertical-align:middle;border:1px solid}
.st.fixed{color:#1b5e20;background:#e8f5e9;border-color:#a5d6a7}.st.open{color:#b71c1c;background:#ffebee;border-color:#ef9a9a}
.st.review{color:#6d4c00;background:#fff8e1;border-color:#ffe082}
.f.fixed h4,.f.fixed pre,.f.fixed p{opacity:.55}.f.fixed h4 b,.f.fixed h4{text-decoration-color:#2e7d32}
details.allfixed>summary{opacity:.6}.tick{color:#2e7d32;font-weight:700}
.prog{margin:1em 0 .5em}.bar{display:flex;height:14px;border-radius:7px;overflow:hidden;background:#eee}
.bar span{display:block}.bar .fixed{background:#43a047}.bar .review{background:#ffca28}.bar .open{background:#e53935}
.sevtab{border-collapse:collapse;margin:.6em 0 1em}.sevtab td,.sevtab th{padding:.15em .7em;border-bottom:1px solid #eee;text-align:right}
.sevtab th:first-child,.sevtab td:first-child{text-align:left}
.filters{position:sticky;top:0;background:#fff;padding:.5em 0;z-index:2;border-bottom:1px solid #eee;margin-bottom:1em}
.filters button{font:inherit;padding:.2em .8em;margin-right:.3em;border:1px solid #ccc;border-radius:6px;background:#fafafa;cursor:pointer}
.filters button.on{background:#1d1d1f;color:#fff;border-color:#1d1d1f}.filters input{font:inherit;padding:.2em .5em;margin-left:.6em;width:14em}
body.show-open .f:not(.open),body.show-fixed .f:not(.fixed),body.show-review .f:not(.review){display:none}
details.hide,h2.hide{display:none}
"""

HTML_SCRIPT = """<script>
(function(){
  var mode='all', q='';
  function apply(){
    document.body.className = mode==='all' ? '' : 'show-'+mode;
    var ds=document.querySelectorAll('details');
    for(var i=0;i<ds.length;i++){
      var d=ds[i], name=(d.getAttribute('data-b')||'').toLowerCase();
      var vis = mode==='all' || d.querySelector('.f.'+mode);
      if(q && name.indexOf(q)<0 && (d.getAttribute('data-ids')||'').indexOf(q)<0) vis=false;
      d.classList.toggle('hide', !vis);
    }
    var hs=document.querySelectorAll('h2[data-cat]');
    for(var j=0;j<hs.length;j++){
      var h=hs[j], n=h.nextElementSibling, any=false;
      while(n && n.tagName==='DETAILS'){ if(!n.classList.contains('hide')) any=true; n=n.nextElementSibling; }
      h.classList.toggle('hide', !any);
    }
  }
  var bs=document.querySelectorAll('.filters button');
  for(var k=0;k<bs.length;k++) bs[k].addEventListener('click', function(){
    for(var m=0;m<bs.length;m++) bs[m].classList.remove('on');
    this.classList.add('on'); mode=this.getAttribute('data-m'); apply();
  });
  document.getElementById('flt').addEventListener('input', function(){ q=this.value.trim().toLowerCase(); apply(); });
})();
</script>"""


def _bin_version():
    """$VersionNumber of the tree the status was taken against."""
    try:
        txt = open(os.path.join(os.path.dirname(BIN), "src", "version.h")).read()
        return re.search(r'MATHILDA_VERSION_STRING "([^"]+)"', txt).group(1)
    except (OSError, AttributeError):
        return "?"


def html_report(rows=None, status=None, adj=None):
    import html as H
    import shutil
    from collections import Counter
    rows = rows or json.load(open(BACKLOG))
    status = status if status is not None else json.load(open(STATUS_JSON))
    if adj is None:
        adj = json.load(open(ADJUDICATION)) if os.path.exists(ADJUDICATION) else {}
    if not os.path.exists(HTML_SOURCE):          # first run: keep the original
        shutil.copyfile(HTML_OUT, HTML_SOURCE)
    src = open(HTML_SOURCE).read()
    byid = {r["id"]: r for r in rows}

    def eff(rid):
        s = effective_status(rid, status, adj)
        return s, HTML_GROUP.get(s, "review")

    tally = Counter()
    for r in rows:
        tally[(r["severity"], eff(r["id"])[1])] += 1
        tally[("all", eff(r["id"])[1])] += 1

    head, *blocks = re.split(r"(?=<details>)", src)
    out_blocks = []
    for b in blocks:
        name = H.unescape(re.search(r"<summary>.*?<b>(.*?)</b>", b, re.S).group(1))
        k = 0
        groups = []

        def stamp(m):
            nonlocal k
            k += 1
            rid = f"{name}#{k}"
            s, g = eff(rid)
            groups.append(g)
            why = ""
            if rid in adj and s == adj[rid]["verdict"]:
                why = adj[rid].get("reason", "")
            note = status.get(rid, {}).get("note", "")
            tip = H.escape(" -- ".join(x for x in (why, note) if x), quote=True)
            label = s.replace("SAME", "OPEN").replace("STILL-BROKEN", "OPEN")
            pill = (f"<span class='st {g}' title='{tip}'>"
                    f"{'&#10003; ' if g == 'fixed' else ''}{label}</span>")
            return (f"<div class='f {g}' id='{H.escape(rid)}'><h4>"
                    f"<small><code>{H.escape(rid)}</code></small> {pill} ")
        b = re.sub(r"<div class='f'><h4>", stamp, b)
        nfix = groups.count("fixed")
        if not groups:
            out_blocks.append(b)
            continue
        cls = " class='allfixed'" if nfix == len(groups) else ""
        summ = (f" <span class='st {'fixed' if nfix == len(groups) else 'open' if nfix == 0 else 'review'}'>"
                f"{'<span class=tick>&#10003;</span> ' if nfix == len(groups) else ''}"
                f"{nfix}/{len(groups)} fixed</span>")
        ids = " ".join(f"{name}#{i}".lower() for i in range(1, len(groups) + 1))
        b = b.replace("<details>", f"<details{cls} data-b='{H.escape(name, quote=True)}' "
                                   f"data-ids='{H.escape(ids, quote=True)}'>", 1)
        b = b.replace("</summary>", summ + "</summary>", 1)
        out_blocks.append(b)

    # category headings get their own fixed / total
    body = "".join(out_blocks)
    head_html, _, first_h2 = head.rpartition("<h2>")
    if not head_html:                            # no heading split: keep as is
        head_html, first_h2 = head, ""
    full = ("<h2>" + first_h2 if first_h2 else "") + body

    def cat_tag(m):
        return f"<h2 data-cat='{H.escape(m.group(1), quote=True)}'>{m.group(1)}</h2>"
    full = re.sub(r"<h2>(.*?)</h2>", cat_tag, full)

    meta = status.get("_meta")
    if meta:
        stamp_txt = f"{meta['time']} against Mathilda v{meta['version']}"
        if meta.get("scope") == "partial":
            stamp_txt += " (last run was partial; other rows keep their earlier status)"
    else:
        mt = time.strftime("%Y-%m-%d %H:%M",
                           time.localtime(os.path.getmtime(STATUS_JSON)))
        stamp_txt = f"{mt} (binary version not recorded)"
    tot = len(rows)
    fx, rv, op = (tally[("all", g)] for g in ("fixed", "review", "open"))
    pct = lambda n: f"{100.0 * n / tot:.2f}%"
    sev_rows = "".join(
        f"<tr><td><span class='sev {sv}'>{sv}</span></td>"
        f"<td>{tally[(sv, 'fixed')]}</td><td>{tally[(sv, 'review')]}</td>"
        f"<td>{tally[(sv, 'open')]}</td>"
        f"<td>{sum(tally[(sv, g)] for g in ('fixed', 'review', 'open'))}</td></tr>"
        for sv in SEV_RANK)
    progress = (
        f"<h2>Fix progress</h2><p>Findings last re-run {stamp_txt}; page rendered "
        f"{time.strftime('%Y-%m-%d %H:%M')} by <code>tools/audit_backlog.py</code>. "
        f"<b>{fx}</b> fixed, <b>{rv}</b> need a verdict (changed since the audit, "
        f"not yet confirmed), <b>{op}</b> still open.</p>"
        f"<div class='prog'><div class='bar'>"
        f"<span class='fixed' style='width:{pct(fx)}' title='fixed {fx}'></span>"
        f"<span class='review' style='width:{pct(rv)}' title='needs verdict {rv}'></span>"
        f"<span class='open' style='width:{pct(op)}' title='open {op}'></span></div></div>"
        f"<table class='sevtab'><tr><th>severity</th><th>fixed</th><th>needs verdict</th>"
        f"<th>open</th><th>total</th></tr>{sev_rows}</table>"
        f"<div class='filters'>"
        f"<button class='on' data-m='all'>All</button>"
        f"<button data-m='open'>Open ({op})</button>"
        f"<button data-m='review'>Needs verdict ({rv})</button>"
        f"<button data-m='fixed'>Fixed ({fx})</button>"
        f"<input id='flt' placeholder='filter: builtin or id'></div>")
    head_html = head_html.replace("</style>", HTML_STYLE + "</style>", 1)
    page = head_html + progress + full
    page = page.replace("</html>", HTML_SCRIPT + "</html>") if "</html>" in page \
        else page + HTML_SCRIPT
    tmp = HTML_OUT + ".tmp"
    open(tmp, "w").write(page)
    os.replace(tmp, HTML_OUT)
    print(f"report.html: {fx} fixed, {rv} need a verdict, {op} open (of {tot})")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "parse":
        rs = parse()
        print(f"{len(rs)} findings -> {BACKLOG}")
    elif cmd == "run":
        run(sys.argv[2:])
    elif cmd == "report":
        report()
    elif cmd == "html":
        html_report()
    elif cmd == "pinned":
        adjudicate_pinned(sys.argv[2:])
    else:
        sys.exit(__doc__)
