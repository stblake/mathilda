#!/usr/bin/env python3
"""Run FeynCalc's upstream unit tests under Mathilda and ratchet the pass set.

    python3 tools/check_feyncalc.py                 # all unit-test directories
    python3 tools/check_feyncalc.py --dir Lorentz --dir Dirac
    python3 tools/check_feyncalc.py --file Contract.test
    python3 tools/check_feyncalc.py --integration   # the IntegrationTests drivers too
    python3 tools/check_feyncalc.py --update-baseline
    python3 tools/check_feyncalc.py --triage        # message / verdict histograms

The vendored package lives in packages/feyncalc/ (local, untracked); without it
this check prints a skip and exits 0.  Each upstream driver
(Tests/<Dir>/<Dir>.mt) is run UNMODIFIED by tools/feyncalc/run_fc_tests.m, a
minimal MUnit`Test.  The driver loads every .test file of its directory and
keeps the groups that match $OnlySubTest (the upstream onlySubTest hook), so a
run is sharded one process per .test file: shards run in parallel, a crash
loses one file's tests rather than a directory's, and FeynCalc loads in about a
second, which is what a shard costs on top of its tests.

Verdicts, per test id: PASS, FAIL (wrong value), MSG (right value, unexpected
messages), ABORT, TIMEOUT, SKIP (typesetting tests: box output is out of
scope), NOTRUN (an id in a .test file that no line reported: its shard crashed
or timed out; an id in a group the upstream driver never runs is not counted).
Results are written to packages/feyncalc/harness/results/
<Dir>.tsv.  The ratchet: every id in harness/BASELINE_PASS.txt that was run
must still PASS (exit 1 otherwise); newly passing ids are listed so the
baseline can be raised with --update-baseline."""
import argparse
import collections
import concurrent.futures as cf
import os
import re
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(REPO, "packages", "feyncalc")
TESTS = os.path.join(PKG, "Tests")
HARNESS = os.path.join(PKG, "harness")
BASELINE = os.path.join(HARNESS, "BASELINE_PASS.txt")
RUNNER = os.path.join(REPO, "tools", "feyncalc", "run_fc_tests.m")
UNIT_DIRS = ["Shared", "NonCommAlgebra", "Lorentz", "Dirac", "SUN", "Pauli",
             "Feynman", "Tables", "LoopIntegrals", "ExportImport"]
VERDICTS = ["PASS", "FAIL", "MSG", "ABORT", "TIMEOUT", "SKIP", "NOTRUN"]

GROUP_RE = re.compile(r"^(Tests`[A-Za-z]+`[A-Za-z0-9]+)\s*=", re.M)
ID_RE = re.compile(r'\{\s*"([A-Za-z][A-Za-z0-9]*-ID[0-9]+)"\s*,')


def drivers(dirs, integration):
    """(label, .mt path, directory of its .test files) for each driver."""
    out = []
    for d in dirs:
        mt = os.path.join(TESTS, d, d + ".mt")
        if os.path.isfile(mt):
            out.append((d, mt, os.path.join(TESTS, d)))
        if integration:
            idir = os.path.join(TESTS, d, "IntegrationTests")
            for f in sorted(os.listdir(idir)) if os.path.isdir(idir) else []:
                if f.endswith(".mt"):
                    out.append((d + "/IntegrationTests", os.path.join(idir, f), idir))
    return out


def shards(label, mt, tdir, only_files, per_dir=False):
    """One shard per .test file: its groups (for $OnlySubTest) and test ids.
    per_dir: one shard for the whole directory (the reference kernel's run,
    where a process per file would pay its ~30 s load 260 times)."""
    out = []
    for f in sorted(os.listdir(tdir)):
        if not f.endswith(".test") or (only_files and f not in only_files):
            continue
        text = open(os.path.join(tdir, f), encoding="utf-8", errors="replace").read()
        groups = GROUP_RE.findall(text)
        if groups:
            out.append({"label": label, "mt": mt, "file": f, "groups": groups,
                        "ids": ID_RE.findall(text)})
    if per_dir and out:
        out = [{"label": label, "mt": mt, "file": "*", "ids": [i for s in out for i in s["ids"]],
                "groups": [g for s in out for g in s["groups"]]}]
    return out


def run_shard(sh, binary, timeout, shard_timeout):
    # `X$ after a backtick: the runner wraps each as .*re.*, and a bare group
    # name would also select every group it is a prefix of. A whole-directory
    # shard selects everything ("" is the upstream "no filter").
    only = "" if sh["file"] == "*" else "|".join("`" + g.split("`")[-1] + "$" for g in sh["groups"])
    env = dict(os.environ, FC_ROOT=PKG, FC_MT=sh["mt"], FC_ONLYSUB=only,
               FC_TIMEOUT=str(timeout))
    t0 = time.time()
    try:
        p = subprocess.run([binary, "-file", RUNNER], cwd=REPO, env=env,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           timeout=shard_timeout)
        out, rc = p.stdout.decode("utf-8", "replace"), p.returncode
    except subprocess.TimeoutExpired as e:
        out, rc = (e.stdout or b"").decode("utf-8", "replace"), "timeout"
    results, details, ended = {}, {}, False
    for line in out.splitlines():
        f = line.split("\t")
        if f[0] == "FC" and len(f) >= 3:
            key = f[2]
            n = 2
            while key in results:          # an id repeated across groups
                key, n = f"{f[2]}#{n}", n + 1
            results[key] = (f[1], micros_to_secs(f[3]) if len(f) > 3 else "")
        elif f[0] == "FCD" and len(f) >= 2:
            details.setdefault(f[1], f[2:])
        elif line == "FCEND":
            ended = True
    # An id no line reported: lost if the shard died before FCEND, otherwise
    # in a group the driver itself never runs (upstream skips it too: the
    # typesetting groups behind $VersionNumber >= 10, the Abort groups a
    # driver does not list), which is not counted.
    unused = 0
    for i in sh["ids"]:
        if i not in results:
            if ended:
                unused += 1
            else:
                results[i] = ("NOTRUN", "")
    msgs = collections.Counter(m.group(1) for m in
                               re.finditer(r"^([A-Za-z$`]+::[A-Za-z0-9]+)", out, re.M))
    return {"shard": sh, "results": results, "details": details, "ended": ended, "unused": unused,
            "rc": rc, "secs": time.time() - t0, "msgs": msgs,
            "tail": "\n".join(out.splitlines()[-5:])}


def micros_to_secs(field):
    """The driver reports integer microseconds; the tsv carries seconds."""
    try:
        return f"{int(field) / 1e6:.4f}"
    except ValueError:
        return field


def id_files(dirs):
    """(dir, id) -> .test file, for results from a per-dir run (file "*")."""
    out = {}
    for lab, mt, td in drivers(dirs, False):
        for sh in shards(lab, mt, td, set()):
            for i in sh["ids"]:
                out.setdefault((lab, i), sh["file"])
    return out


def read_results(rdir):
    """{(dir, id): (verdict, seconds, file)} from a results directory."""
    out = {}
    for f in sorted(os.listdir(rdir)) if os.path.isdir(rdir) else []:
        if not f.endswith(".tsv"):
            continue
        lab = f[:-4]
        for line in open(os.path.join(rdir, f), encoding="utf-8").read().splitlines()[1:]:
            c = line.split("\t")
            if len(c) >= 4:
                try:
                    secs = float(c[2])
                except ValueError:
                    secs = 0.0
                out[(lab, c[0])] = (c[1], secs, c[3])
    return out


def timing_report(mine_dir, ref_dir, dirs, out_path, floor=1e-3):
    """Speed of this kernel against a reference run (both from the same harness).

    The headline is the geometric mean, over .test files, of T_ref / T_mine,
    counting only ids that PASS in both runs: a wrong answer is not timed.
    Each file's total is floored at `floor` seconds on both sides, so a file
    of trivial tests reads as parity rather than as noise.  > 1 means this
    kernel is faster."""
    import math
    mine, ref = read_results(mine_dir), read_results(ref_dir)
    files = id_files(dirs)
    both = [k for k in mine if k in ref and mine[k][0] == "PASS" and ref[k][0] == "PASS"]
    per_file = collections.defaultdict(lambda: [0.0, 0.0, 0])
    per_dir = collections.defaultdict(lambda: [0.0, 0.0, 0])
    for k in both:
        fn = files.get(k, mine[k][2])
        for agg in (per_file[(k[0], fn)], per_dir[k[0]]):
            agg[0] += mine[k][1]
            agg[1] += ref[k][1]
            agg[2] += 1
    def ratio(m, r):
        return max(r, floor) / max(m, floor)
    lines = []
    w = lines.append
    logs = [math.log(ratio(m, r)) for m, r, _ in per_file.values()]
    gm = math.exp(sum(logs) / len(logs)) if logs else float("nan")
    tm = sum(v[0] for v in per_dir.values())
    tr = sum(v[1] for v in per_dir.values())
    def show(d):
        r = os.path.relpath(d, REPO)
        return d if r.startswith("..") else r
    w(f"# FeynCalc timing: {show(mine_dir)} vs {show(ref_dir)}")
    w("")
    w(f"Speed-up = T_ref / T_mine (> 1: this kernel is faster); ids passing in both: {len(both)}, "
      f"files: {len(per_file)}.")
    w("")
    w(f"- **geomean over files: {gm:.3f}x**")
    w(f"- total compute: mine {tm:.1f} s, ref {tr:.1f} s, ratio {ratio(tm, tr):.3f}x")
    w("")
    w("| directory | ids | mine s | ref s | speed-up |")
    w("|---|---:|---:|---:|---:|")
    for d in sorted(per_dir):
        m, r, n = per_dir[d]
        w(f"| {d} | {n} | {m:.2f} | {r:.2f} | {ratio(m, r):.3f}x |")
    rows = [(k, v) for k, v in per_file.items()]
    w("")
    w("## Worst 30 files by speed-up")
    w("")
    w("| file | ids | mine s | ref s | speed-up |")
    w("|---|---:|---:|---:|---:|")
    for (d, fn), (m, r, n) in sorted(rows, key=lambda kv: ratio(kv[1][0], kv[1][1]))[:30]:
        w(f"| {d}/{fn} | {n} | {m:.3f} | {r:.3f} | {ratio(m, r):.3f}x |")
    w("")
    w("## Top 30 files by this kernel's time")
    w("")
    w("| file | ids | mine s | ref s | speed-up |")
    w("|---|---:|---:|---:|---:|")
    for (d, fn), (m, r, n) in sorted(rows, key=lambda kv: -kv[1][0])[:30]:
        w(f"| {d}/{fn} | {n} | {m:.3f} | {r:.3f} | {ratio(m, r):.3f}x |")
    w("")
    w("## Top 30 ids by this kernel's time (any verdict)")
    w("")
    w("| id | verdict | mine s | ref verdict | ref s |")
    w("|---|---|---:|---|---:|")
    for k in sorted(mine, key=lambda k: -mine[k][1])[:30]:
        rv = ref.get(k, ("-", 0.0))
        w(f"| {k[0]}/{k[1]} | {mine[k][0]} | {mine[k][1]:.3f} | {rv[0]} | {rv[1]:.3f} |")
    text = "\n".join(lines) + "\n"
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(text)
    print("\n".join(lines[:12 + len(per_dir)]))
    print(f"timing report written to {os.path.relpath(out_path, REPO)}")


def load_baseline():
    base = set()
    if os.path.isfile(BASELINE):
        for line in open(BASELINE, encoding="utf-8"):
            line = line.rstrip("\n")
            if line and not line.startswith("#"):
                base.add(tuple(line.split("\t", 1)))
    return base


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dir", action="append", help="test directory (repeatable)")
    ap.add_argument("--file", action="append", help="only this .test file (repeatable)")
    ap.add_argument("--integration", action="store_true", help="include IntegrationTests drivers")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--timeout", type=int, default=60, help="per-test budget, seconds")
    ap.add_argument("--shard-timeout", type=int, default=3600, help="per-process budget, seconds")
    ap.add_argument("--binary", default=os.path.join(REPO, "Mathilda"))
    ap.add_argument("--per-dir", action="store_true",
                    help="one process per directory (for the reference kernel: --binary wolframscript)")
    ap.add_argument("--no-ratchet", action="store_true", help="report only; never fail on the baseline")
    ap.add_argument("--results-dir", default=os.path.join(HARNESS, "results"),
                    help="where the per-directory .tsv go (an A/B run against another binary)")
    ap.add_argument("--update-baseline", action="store_true")
    ap.add_argument("--triage", action="store_true", help="print verdict and message histograms")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--compare-timing", metavar="REF_DIR",
                    help="after the run, compare per-test times against a reference results "
                         "directory (e.g. a --per-dir wolframscript run); writes harness/TIMING.md")
    ap.add_argument("--report-only", action="store_true",
                    help="with --compare-timing: skip the run, compare the existing --results-dir")
    a = ap.parse_args()

    if not os.path.isdir(TESTS):
        print("check-feyncalc: packages/feyncalc not present -- skipped")
        return 0
    dirs = a.dir or UNIT_DIRS
    if a.report_only:
        if not a.compare_timing:
            print("check-feyncalc: --report-only needs --compare-timing")
            return 2
        timing_report(a.results_dir, a.compare_timing, dirs, os.path.join(HARNESS, "TIMING.md"))
        return 0
    work = [s for (lab, mt, td) in drivers(dirs, a.integration)
            for s in shards(lab, mt, td, set(a.file or []), a.per_dir)]
    if not work:
        print("check-feyncalc: nothing selected")
        return 2
    # Longest first: the slowest shards should not start last.
    work.sort(key=lambda s: -len(s["ids"]))
    t0 = time.time()
    runs = []
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(run_shard, s, a.binary, a.timeout, a.shard_timeout) for s in work]
        for fu in cf.as_completed(futs):
            r = fu.result()
            runs.append(r)
            if not a.quiet:
                c = collections.Counter(v for v, _ in r["results"].values())
                bad = "" if r["ended"] else f"  [no FCEND, rc={r['rc']}]"
                print(f"  {r['shard']['label']}/{r['shard']['file']}: "
                      f"{c['PASS']}/{len(r['results'])} in {r['secs']:.0f}s{bad}", flush=True)

    # Results per label.
    os.makedirs(a.results_dir, exist_ok=True)
    per = collections.defaultdict(dict)
    for r in runs:
        for i, (v, secs) in r["results"].items():
            d = r["details"].get(i, ["", "", ""])
            per[r["shard"]["label"]][i] = (v, secs, r["shard"]["file"], *d[:3])
    for lab, res in per.items():
        path = os.path.join(a.results_dir, lab.replace("/", "_") + ".tsv")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("id\tverdict\tseconds\tfile\tgot\texpected\tmessages\n")
            for i in sorted(res):
                fh.write(i + "\t" + "\t".join(str(x) for x in res[i]) + "\n")

    # Summary.
    print(f"\n{'directory':<28}" + "".join(f"{v:>8}" for v in VERDICTS) + f"{'total':>8}{'pass%':>8}")
    tot = collections.Counter()
    for lab in sorted(per):
        c = collections.Counter(v[0] for v in per[lab].values())
        tot.update(c)
        n = sum(c.values())
        print(f"{lab:<28}" + "".join(f"{c[v]:>8}" for v in VERDICTS) +
              f"{n:>8}{100.0 * c['PASS'] / max(1, n - c['SKIP']):>7.1f}%")
    n = sum(tot.values())
    print(f"{'ALL':<28}" + "".join(f"{tot[v]:>8}" for v in VERDICTS) +
          f"{n:>8}{100.0 * tot['PASS'] / max(1, n - tot['SKIP']):>7.1f}%")
    print(f"wall {time.time() - t0:.0f}s, {len(runs)} shards, {a.jobs} jobs; "
          f"{sum(r['unused'] for r in runs)} ids in groups the drivers do not run")

    if a.triage:
        msgs = collections.Counter()
        for r in runs:
            msgs.update(r["msgs"])
        print("\nmessages printed (all shards):")
        for m, k in msgs.most_common(40):
            print(f"  {k:>6}  {m}")
        crashed = [r for r in runs if not r["ended"]]
        if crashed:
            print("\nshards without FCEND:")
            for r in crashed:
                print(f"  {r['shard']['label']}/{r['shard']['file']} rc={r['rc']}: {r['tail'][-300:]}")

    if a.compare_timing:
        print()
        timing_report(a.results_dir, a.compare_timing, dirs, os.path.join(HARNESS, "TIMING.md"))

    # Ratchet.
    now_pass = {(lab, i) for lab, res in per.items() for i, v in res.items() if v[0] == "PASS"}
    ran = {(lab, i) for lab, res in per.items() for i in res}
    base = load_baseline()
    regressed = sorted(k for k in base if k in ran and k not in now_pass)
    gained = sorted(now_pass - base)
    if gained:
        print(f"\n{len(gained)} newly passing (not in the baseline)" +
              ("" if a.update_baseline else "; raise it with --update-baseline"))
    if a.update_baseline:
        keep = {k for k in base if k not in ran}
        os.makedirs(HARNESS, exist_ok=True)
        with open(BASELINE, "w", encoding="utf-8") as fh:
            fh.write("# FeynCalc unit tests that pass under Mathilda (tools/check_feyncalc.py)\n")
            for lab, i in sorted(keep | now_pass):
                fh.write(f"{lab}\t{i}\n")
        print(f"baseline: {len(keep | now_pass)} passing ids written to {os.path.relpath(BASELINE, REPO)}")
        return 0
    if regressed and not a.no_ratchet:
        print(f"\nREGRESSED: {len(regressed)} baseline ids no longer pass:")
        for lab, i in regressed[:50]:
            v = per[lab].get(i, ("?",))[0]
            print(f"  {lab}\t{i}\t{v}")
        return 1
    print("\nratchet OK" + ("" if base else " (no baseline yet)"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
