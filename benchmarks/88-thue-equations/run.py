#!/usr/bin/env python3
"""benchmarks/88-thue-equations/run.py -- the Thue stress harness.

For every case in cases.py: solve it in Mathilda's Solve[..., Integers] AND in
PARI/GP's thue() (the gold-standard oracle), time both, compare the finite
solution sets, and classify.  Writes REPORT.md + results.json.

Verdicts
    CORRECT     both solved, sets identical
    WRONG       both solved, sets DIFFER            <- COMPLETENESS BUG
    CRASH       Mathilda errored / no output        <- BUG
    TIMEOUT     Mathilda exceeded the time budget   <- bottleneck
    DECLINE     Mathilda unevaluated, PARI solved   <- honest coverage gap
    UNVERIFIED  Mathilda solved, PARI could not cross-check

Any WRONG or CRASH makes the run exit nonzero.

    python3 run.py                 # all cases
    python3 run.py --only binom    # label substring filter
    MATHILDA_BIN=../../Mathilda MATH_TIMEOUT=30 PARI_TIMEOUT=25 python3 run.py
"""

import ast
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cases import CASES  # noqa: E402

REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
MATHILDA = os.environ.get("MATHILDA_BIN", os.path.join(REPO, "Mathilda"))
GP = os.environ.get("GP_BIN", "gp")
MATH_TIMEOUT = float(os.environ.get("MATH_TIMEOUT", "45"))
PARI_TIMEOUT = float(os.environ.get("PARI_TIMEOUT", "30"))


# ----------------------------------------------------------------------------
#  Mathilda side
# ----------------------------------------------------------------------------
def run_mathilda(case):
    form, m = case["form"], case["m"]
    script = (
        f"tr = Timing[Solve[{form} == {m} && Element[{{x, y}}, Integers], {{x, y}}, Integers]];\n"
        f"res = Last[tr];\n"
        f'Print["TIME:", First[tr]];\n'
        f'Print["STATUS:", If[Head[res] === List, "SOLVED", "DECLINE"]];\n'
        f'If[Head[res] === List, Print["SOLS:", Sort[Map[({{x, y}} /. #) &, res]]]];\n'
    )
    tmp = os.path.join(HERE, "_case.m")
    with open(tmp, "w") as f:
        f.write(script)
    t0 = time.perf_counter()
    try:
        p = subprocess.run([MATHILDA, tmp], capture_output=True, text=True,
                           timeout=MATH_TIMEOUT)
    except subprocess.TimeoutExpired:
        return dict(status="TIMEOUT", sols=None, ctime=None, wall=MATH_TIMEOUT)
    wall = time.perf_counter() - t0
    out = p.stdout
    if p.returncode != 0 or "STATUS:" not in out:
        return dict(status="CRASH", sols=None, ctime=None, wall=wall,
                    err=(p.stderr or out)[:400])
    ctime = None
    mt = re.search(r"TIME:\s*([0-9.eE+-]+)", out)
    if mt:
        try: ctime = float(mt.group(1))
        except ValueError: ctime = None
    if "STATUS:DECLINE" in out or "STATUS: DECLINE" in out:
        return dict(status="DECLINE", sols=None, ctime=ctime, wall=wall)
    ms = re.search(r"SOLS:\s*(\{.*\})", out, re.S)
    if not ms:
        return dict(status="CRASH", sols=None, ctime=ctime, wall=wall, err="no SOLS line")
    sols = parse_braces(ms.group(1))
    return dict(status="SOLVED", sols=sols, ctime=ctime, wall=wall)


def parse_braces(s):
    """'{{-1, -1}, {1, 0}}' or '{}' -> sorted set of (x,y) int tuples, or None."""
    s = s.strip().replace("{", "[").replace("}", "]")
    try:
        v = ast.literal_eval(s)
    except (ValueError, SyntaxError):
        return None
    if not isinstance(v, list):
        return None
    out = set()
    for t in v:
        if isinstance(t, (list, tuple)) and len(t) == 2 and all(isinstance(z, int) for z in t):
            out.add((t[0], t[1]))
        else:
            return None
    return frozenset(out)


# ----------------------------------------------------------------------------
#  PARI oracle
# ----------------------------------------------------------------------------
def run_pari(case):
    P, m = case["P"], case["m"]
    # gettime() around thueinit+thue = PARI COMPUTE time (excludes gp's ~30ms
    # process startup), so the ratio vs Mathilda's Timing[] is apples-to-apples.
    cmd = (f"gettime();r = thue(thueinit({P}, 0), {m});t=gettime();"
           f'print("PARITIME:", t);print("PARISOLS:", r);')
    t0 = time.perf_counter()
    try:
        p = subprocess.run([GP, "-q", "-s", "800000000"], input=cmd,
                           capture_output=True, text=True, timeout=PARI_TIMEOUT)
    except subprocess.TimeoutExpired:
        return dict(status="TIMEOUT", sols=None, wall=PARI_TIMEOUT, ctime=None)
    except FileNotFoundError:
        return dict(status="NOGP", sols=None, wall=0.0, ctime=None)
    wall = time.perf_counter() - t0
    ct = re.search(r"PARITIME:\s*([0-9]+)", p.stdout)
    ctime = (int(ct.group(1)) / 1000.0) if ct else None       # gettime() is ms
    mo = re.search(r"PARISOLS:\s*(\[.*\])", p.stdout, re.S)
    if not mo or p.returncode != 0:
        return dict(status="ERROR", sols=None, wall=wall, ctime=ctime, err=(p.stderr or p.stdout)[:300])
    try:
        v = ast.literal_eval(mo.group(1).replace("\n", ""))
        sols = frozenset((int(a), int(b)) for a, b in v)
    except Exception:
        return dict(status="ERROR", sols=None, wall=wall, ctime=ctime, err="parse")
    return dict(status="OK", sols=sols, wall=wall, ctime=ctime)


# ----------------------------------------------------------------------------
#  classify + report
# ----------------------------------------------------------------------------
def classify(mo, po):
    ms, ps = mo["status"], po["status"]
    if ms == "TIMEOUT": return "TIMEOUT"
    if ms == "CRASH":   return "CRASH"
    if ms == "SOLVED":
        if ps == "OK":
            return "CORRECT" if mo["sols"] == po["sols"] else "WRONG"
        return "UNVERIFIED"
    # DECLINE
    return "DECLINE"


def main():
    only = None
    for a in sys.argv[1:]:
        if a.startswith("--only"):
            only = a.split("=", 1)[1] if "=" in a else sys.argv[sys.argv.index(a) + 1]
    cases = [c for c in CASES if (only is None or only in c["label"])]

    results = []
    counts = {}
    print(f"Running {len(cases)} Thue cases (Mathilda vs PARI oracle)...\n")
    for i, c in enumerate(cases, 1):
        mo = run_mathilda(c)
        po = run_pari(c)
        verdict = classify(mo, po)
        counts[verdict] = counts.get(verdict, 0) + 1
        mn = len(mo["sols"]) if mo.get("sols") is not None else "-"
        pn = len(po["sols"]) if po.get("sols") is not None else "-"
        mt = f'{mo["ctime"]*1000:.1f}ms' if mo.get("ctime") else (f'{mo["wall"]:.1f}s*' if mo.get("wall") else "-")
        flag = "  <<<" if verdict in ("WRONG", "CRASH") else ""
        print(f"[{i:3d}/{len(cases)}] {verdict:10s} {c['label']:28s} "
              f"M:{mn} P:{pn}  {mt}{flag}")
        if verdict == "WRONG":
            print(f"        Mathilda: {sorted(mo['sols'])}")
            print(f"        PARI    : {sorted(po['sols'])}")
        if verdict == "CRASH" and mo.get("err"):
            print(f"        err: {mo['err'][:200]}")
        results.append(dict(label=c["label"], family=c["family"], note=c["note"],
                            form=c["form"], m=c["m"], verdict=verdict,
                            m_status=mo["status"], m_ctime=mo.get("ctime"),
                            m_wall=mo.get("wall"), m_count=(len(mo["sols"]) if mo.get("sols") is not None else None),
                            p_status=po["status"], p_wall=po.get("wall"), p_ctime=po.get("ctime"),
                            p_count=(len(po["sols"]) if po.get("sols") is not None else None),
                            m_err=mo.get("err")))

    write_report(results, counts)
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(dict(counts=counts, results=results), f, indent=1)
    try:
        os.remove(os.path.join(HERE, "_case.m"))   # scratch file the runner writes
    except OSError:
        pass

    print("\n" + "  ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    bugs = counts.get("WRONG", 0) + counts.get("CRASH", 0)
    if bugs:
        print(f"\n!!! {bugs} BUG(S) FOUND (WRONG/CRASH) -- see REPORT.md")
    return 1 if bugs else 0


def write_report(results, counts):
    lines = []
    L = lines.append
    L("# Benchmark 88 — Thue-equation stress tests\n")
    L("Mathilda's `Solve[F(x,y)==m && Element[{x,y},Integers], {x,y}, Integers]` "
      "vs **PARI/GP `thue()`** (gold-standard oracle) over an adversarial corpus. "
      "A `WRONG` verdict (finite sets differ) is a completeness bug.\n")
    order = ["CORRECT", "WRONG", "CRASH", "TIMEOUT", "DECLINE", "UNVERIFIED"]
    L("| verdict | n | meaning |")
    L("|---|---|---|")
    meaning = dict(CORRECT="matches PARI", WRONG="**differs from PARI — BUG**",
                   CRASH="**errored / no output — BUG**", TIMEOUT="exceeded budget — bottleneck",
                   DECLINE="unevaluated (honest gap); PARI solved", UNVERIFIED="PARI could not cross-check")
    for k in order:
        if counts.get(k):
            L(f"| {k} | {counts[k]} | {meaning[k]} |")
    L("")

    bugs = [r for r in results if r["verdict"] in ("WRONG", "CRASH")]
    if bugs:
        L("## 🐞 BUGS\n")
        for r in bugs:
            L(f"- **{r['verdict']}** `{r['label']}` — `{r['form']} == {r['m']}` "
              f"(Mathilda {r['m_count']} vs PARI {r['p_count']}). {r['note']}")
        L("")

    solved = [r for r in results if r["verdict"] == "CORRECT" and r["m_ctime"]]
    solved.sort(key=lambda r: -(r["m_ctime"] or 0))
    if solved:
        L("## Bottlenecks — slowest CORRECT solves\n")
        L("Both times are COMPUTE time (Mathilda `Timing[]`, PARI `gettime()`), "
          "excluding each process's ~30–50 ms startup, so the ratio is fair.\n")
        L("| label | form == m | Mathilda | PARI | ratio | #sol |")
        L("|---|---|---:|---:|---:|---:|")
        for r in solved[:15]:
            mt = (r["m_ctime"] or 0) * 1000
            pt = (r.get("p_ctime") or 0) * 1000
            ratio = f"{mt/pt:.1f}x" if pt > 0.5 else "-"
            L(f"| {r['label']} | `{r['form']} == {r['m']}` | {mt:.1f}ms | {pt:.0f}ms | {ratio} | {r['m_count']} |")
        L("")

    # per-family coverage
    fams = {}
    for r in results:
        f = r["family"]; fams.setdefault(f, {}).setdefault(r["verdict"], 0)
        fams[f][r["verdict"]] += 1
    L("## Coverage by family\n")
    L("| family | CORRECT | DECLINE | WRONG | CRASH | TIMEOUT |")
    L("|---|---:|---:|---:|---:|---:|")
    for f in sorted(fams):
        d = fams[f]
        L(f"| {f} | {d.get('CORRECT',0)} | {d.get('DECLINE',0)} | "
          f"{d.get('WRONG',0)} | {d.get('CRASH',0)} | {d.get('TIMEOUT',0)} |")
    L("")

    # full table
    L("## All cases\n")
    L("| label | form == m | verdict | M#/P# | M time | note |")
    L("|---|---|---|---|---:|---|")
    for r in results:
        mt = f'{r["m_ctime"]*1000:.1f}ms' if r["m_ctime"] else "-"
        L(f"| {r['label']} | `{r['form']} == {r['m']}` | {r['verdict']} | "
          f"{r['m_count']}/{r['p_count']} | {mt} | {r['note']} |")
    with open(os.path.join(HERE, "REPORT.md"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
