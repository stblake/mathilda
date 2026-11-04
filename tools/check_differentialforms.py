#!/usr/bin/env python3
"""check_differentialforms.py -- run the DifferentialForms package corpus.

Issue #125: the DifferentialForms package (Zizza/Jentschura, v3.1) is a
third-party package that exercises package loading, context resolution,
UpValues, legacy syntax (Literal, positional ReplacePart), Dt, Integrate and
Simplify. The package itself is copyrighted and lives untracked in
packages/differentialforms/ (like packages/feyncalc/); the corpus
(tests/differentialforms/corpus.m) and the driver
(tools/differentialforms/run_df_cases.m) are ours.

Usage
    python3 tools/check_differentialforms.py --record --binary B  # record expected answers from kernel B
    python3 tools/check_differentialforms.py               # run Mathilda, compare, ratchet
    python3 tools/check_differentialforms.py --only id1,id2 --verbose
    python3 tools/check_differentialforms.py --update-baseline

Every case is one of
    recorded  compared, semantically, with the recorded expected answer
    identity  must evaluate to True (a theorem the package must satisfy)
    expect    compared with an expected value written in the corpus
and is reported PASS, FAIL, MSG (right answer, unexpected message), TIMEOUT,
NOREF (no recorded answer) or CRASH (the process died on it).

The run is one process; if it dies, the cases it did not reach are rerun one
process each so a single crash cannot hide the rest of the corpus.

Exit status: 0 when no baseline PASS regressed (or the package is absent, in
which case the check is skipped), 1 otherwise.
"""

import argparse
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG_DIR = os.path.join(ROOT, "packages", "differentialforms")
PKG = os.path.join(PKG_DIR, "DifferentialForms.m")
HARNESS = os.path.join(PKG_DIR, "harness")
CORPUS = os.path.join(ROOT, "tests", "differentialforms", "corpus.m")
DRIVER = os.path.join(ROOT, "tools", "differentialforms", "run_df_cases.m")


def run_driver(binary, env_extra, timeout):
    """Run the driver under `binary`; return (stdout lines, finished?)."""
    env = dict(os.environ)
    env.update({"DF_PKG": PKG, "DF_CORPUS": CORPUS})
    env.update(env_extra)
    try:
        p = subprocess.run([binary, "-file", DRIVER], env=env, cwd=ROOT,
                           capture_output=True, text=True, timeout=timeout,
                           errors="replace")
        out = p.stdout
    except subprocess.TimeoutExpired as e:
        out = e.stdout or ""
        if isinstance(out, bytes):
            out = out.decode("utf-8", "replace")
    lines = out.splitlines()
    return lines, any(l.startswith("DFEND") for l in lines)


def corpus_ids():
    """Case ids in corpus order (parsed from the corpus text)."""
    ids = []
    with open(CORPUS, encoding="utf-8") as f:
        for line in f:
            s = line.lstrip()
            if s.startswith('{"'):
                ids.append(s[2:s.index('"', 2)])
    return ids


def m_string(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def record(args):
    lines, done = run_driver(args.binary, {"DF_TIMEOUT": str(args.timeout)},
                             args.timeout * 400)
    rows = [l.split("\t") for l in lines if l.startswith("DFR\t")]
    if not done:
        print("recording run did not finish; %d cases dumped" % len(rows))
    os.makedirs(HARNESS, exist_ok=True)
    with open(os.path.join(HARNESS, "expected.tsv"), "w", encoding="utf-8") as f:
        for r in rows:
            f.write("\t".join(r[1:]) + "\n")
    with open(os.path.join(HARNESS, "expected.m"), "w", encoding="utf-8") as f:
        f.write("$DFExpected = <|\n")
        ents = []
        for r in rows:
            _, cid, _kind, full, msgs = (r + [""] * 5)[:5]
            ml = ", ".join(m_string(m) for m in msgs.split(";") if m)
            ents.append("  %s -> {%s, {%s}}" % (m_string(cid), m_string(full), ml))
        f.write(",\n".join(ents) + "\n|>;\n")
    print("recorded %d expected answers in %s" % (len(rows), HARNESS))
    return 0


def parse_df(lines):
    res = {}
    for l in lines:
        if l.startswith("DF\t"):
            parts = (l.split("\t") + [""] * 6)[:6]
            _, verdict, cid, kind, got, msgs = parts
            res[cid] = (verdict, kind, got, msgs)
    return res


def compare(args):
    exp_m = os.path.join(HARNESS, "expected.m")
    if not os.path.exists(exp_m):
        print("no recorded answers (%s); run with --record --binary B first" % exp_m)
        return 1
    ids = corpus_ids()
    if args.only:
        ids = [i for i in ids if i in set(args.only.split(","))]
    env = {"DF_TIMEOUT": str(args.timeout), "DF_EXPECTED": exp_m}
    if args.only:
        env["DF_ONLY"] = ",".join(ids)
    lines, done = run_driver(args.binary, env, args.timeout * max(len(ids), 1) + 120)
    res = parse_df(lines)
    if not done:
        missing = [i for i in ids if i not in res]
        print("run died after %d cases; rerunning %d one by one" % (len(res), len(missing)))
        for cid in missing:
            e1 = dict(env, DF_ONLY=cid)
            l1, _ = run_driver(args.binary, e1, args.timeout + 120)
            r1 = parse_df(l1)
            res[cid] = r1.get(cid, ("CRASH", "?", "", ""))

    counts = {}
    os.makedirs(HARNESS, exist_ok=True)
    with open(os.path.join(HARNESS, "results.tsv"), "w", encoding="utf-8") as f:
        for cid in ids:
            verdict, kind, got, msgs = res.get(cid, ("CRASH", "?", "", ""))
            counts[verdict] = counts.get(verdict, 0) + 1
            f.write("\t".join([cid, verdict, kind, got, msgs]) + "\n")
            if verdict != "PASS" and (args.verbose or args.only):
                print("%-8s %-24s %s  [%s]" % (verdict, cid, got[:300], msgs))
    total = len(ids)
    npass = counts.get("PASS", 0)
    print("DifferentialForms corpus: %d/%d PASS  %s" % (
        npass, total, "  ".join("%s %d" % kv for kv in sorted(counts.items()) if kv[0] != "PASS")))

    base_f = os.path.join(HARNESS, "BASELINE_PASS.txt")
    passing = {cid for cid in ids if res.get(cid, ("",))[0] == "PASS"}
    if args.update_baseline:
        with open(base_f, "w") as f:
            f.write("\n".join(sorted(passing)) + "\n")
        print("baseline updated: %d passes" % len(passing))
        return 0
    if args.only or not os.path.exists(base_f):
        return 0
    with open(base_f) as f:
        base = {l.strip() for l in f if l.strip()}
    lost = sorted(i for i in base & set(ids) if i not in passing)
    gained = sorted(passing - base)
    if gained:
        print("newly passing (%d): %s" % (len(gained), " ".join(gained)))
    if lost:
        print("REGRESSED (%d): %s" % (len(lost), " ".join(lost)))
        return 1
    print("ratchet OK (baseline %d)" % len(base))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--record", action="store_true",
                    help="run the corpus under --binary and record its answers as expected")
    ap.add_argument("--binary", default=os.path.join(ROOT, "Mathilda"),
                    help="kernel binary to run (accepts `-file script`)")
    ap.add_argument("--only", default="", help="comma-separated case ids")
    ap.add_argument("--timeout", type=int, default=60, help="per-case seconds")
    ap.add_argument("--verbose", action="store_true", help="print every non-PASS case")
    ap.add_argument("--update-baseline", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(PKG):
        print("check-differentialforms: SKIP (%s not present; it is not vendored "
              "in git -- copy DifferentialForms.m there to run the check)" % PKG)
        return 0
    if args.record and os.path.abspath(args.binary) == os.path.join(ROOT, "Mathilda"):
        print("--record needs --binary pointing at the kernel whose answers are expected")
        return 2
    return record(args) if args.record else compare(args)


if __name__ == "__main__":
    sys.exit(main())
