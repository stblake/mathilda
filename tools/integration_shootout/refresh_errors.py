#!/usr/bin/env python3
"""refresh_errors.py -- re-run the FriCAS and Axiom ERROR records of a benchmark, to record
each one's own error message.

    python3 tools/integration_shootout/refresh_errors.py --bench benchmarks/263-hard-integrals-ii [--jobs 8]

Until 2026-11-25 the adapter (adapters.run_fricas_one) recorded every error that
integrate RAISED as "integrate failed": the message FriCAS printed, typically
"integrate: implementation incomplete (has polynomial part)", was in the session output
and was dropped.  A second fault made some runs fail on the adapter's own output file
("File is not writable").  Both are fixed in the adapter; this re-runs only the affected
records, so nothing else in the benchmark has to be recomputed (a full shootout would
re-run every Mathilda method too, since the binary changed).

For every case whose latest saved run (results/<date>.json) has FriCAS or Axiom at ERROR,
the system is run again through the fixed adapter.  The new record is appended to
results/cache.jsonl under the same key, and the case's row in results/<date>.json is
updated: its raw record, and its classification by shootout.classify_run.  A re-run that
now answers is classified like any answer (a load-dependent failure of the first run);
one that answers but fails its own check is judged in Mathilda by the shootout's
cross-check (shootout.crosscheck), exactly as in a full run.  results/<date>.json records the refresh under
"errors_refreshed".
"""
import argparse
import collections
import datetime
import json
import os
import re
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import adapters as A  # noqa: E402
import shootout as SH  # noqa: E402

SYSTEMS = {"FriCAS": "fricas", "Axiom": "axiom"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=30)
    a = ap.parse_args()
    bench = os.path.abspath(a.bench)
    rdir = os.path.join(bench, "results")
    date = sorted(f[:-5] for f in os.listdir(rdir) if re.match(r"\d{4}-\d\d-\d\d\.json$", f))[-1]
    rpath = os.path.join(rdir, date + ".json")
    saved = json.load(open(rpath))
    cases = {c["id"]: c for c in SH.load_corpus(os.path.join(bench, "corpus.tsv"))}
    work = [(r, lab) for r in saved["rows"] for lab in SYSTEMS
            if r["runs"].get(lab) == "ERROR" and r["id"] in cases]
    print(f"{os.path.basename(bench)}: {len(work)} FriCAS/Axiom ERROR records to re-run (run {date})", flush=True)
    tmpdir = tempfile.mkdtemp(prefix="refresh-")
    cache = SH.Cache(os.path.join(rdir, "cache.jsonl"), False)

    def run(item):
        r, lab = item
        c = cases[r["id"]]
        return item, A.run_fricas_one(c, a.timeout, tmpdir, SYSTEMS[lab])

    def ckey(s, m, c):
        return "|".join([s, m or "", c["id"], c["key"], ""])

    done = []
    with ThreadPoolExecutor(a.jobs) as ex:
        for n, ((r, lab), rec) in enumerate(ex.map(run, work), 1):
            cache.put(ckey(SYSTEMS[lab], None, cases[r["id"]]), rec)
            done.append((r, lab, rec))
            if n % 100 == 0:
                print(f"  {n}/{len(work)}", flush=True)
    # an answer that fails its own check is judged in Mathilda, as in the shootout
    xc = SH.crosscheck([cases[r["id"]] for r, _, _ in done], sorted(set(SYSTEMS[l] for _, l, _ in done)),
                       set(), cache, ckey, SH.mathilda_stamp(), tmpdir, a.jobs)
    after = collections.Counter()
    for r, lab, rec in done:
        xv = xc.get(ckey(SYSTEMS[lab], None, cases[r["id"]]))
        r["runs"][lab] = SH.classify_run(SYSTEMS[lab], rec, xv)
        r["raw"][lab] = {k: rec.get(k) for k in ("status", "time", "nums", "cert", "detail")}
        if xv is not None:
            r["raw"][lab]["pzq"] = xv
        r["raw"][lab]["text"] = (rec.get("text") or "")[:300]
        after[(lab, r["runs"][lab])] += 1
    saved.setdefault("errors_refreshed", []).append(
        {"date": datetime.date.today().isoformat(), "records": len(work),
         "after": {f"{k[0]} {k[1]}": v for k, v in sorted(after.items())}})
    json.dump(saved, open(rpath, "w"), indent=0)
    print("after:", dict(after))


if __name__ == "__main__":
    main()
