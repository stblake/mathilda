#!/usr/bin/env python3
"""pzq_verify.py -- verify every answer of a benchmark run by PossibleZeroQ, and nothing else.

    python3 tools/integration_shootout/pzq_verify.py --bench benchmarks/263-hard-integrals-ii [--jobs 8]

Until 2026-11-26 an answer counted as verified when its own system's check passed:
D[F, x] - f evaluated at three rational points to 40 digits, in Mathilda for Mathilda's
answers, in FriCAS / Axiom by complexNumeric, in Maxima and sympy by themselves.
PossibleZeroQ judged only the answers that failed that check.  A point check is
numeric evidence, not a proof, and it held Mathilda's own answers to a weaker test than
the one that decided "wrong" for everyone else.

Now every answer of every system is judged by the one test, in Mathilda:

    PossibleZeroQ[D[F, x] - f]    (crosscheck.run_batch: 60 s, and a True counts only
                                   where D[F, x] - f evaluates to a number at two
                                   points -- PossibleZeroQ answers True for an
                                   expression it cannot evaluate)

    zero      verified      (ELEM, or SPECIAL when the answer has a special function)
    nonzero   WRONG
    anything else -- undecided, the answer does not parse, does not translate into
              Mathilda (crosscheck.to_wl), or its text was cut when the run was
              recorded -- UNVERIFIED, never verified on numeric evidence

The answers are the full texts of results/cache.jsonl (results/<date>.json keeps 300
characters).  Mathilda's and Wolfram's answers are read as they are; FriCAS, Axiom,
Maxima and sympy answers are translated by crosscheck.to_wl.  Verdicts are cached in
results/pzq.jsonl, keyed by the answer text and the Mathilda build (shootout.mathilda_stamp:
the binary's and the .m packages' mtimes), so a rerun judges only new answers.  The case rows of results/<date>.json are reclassified in place; each
answer's raw record gets "pzq" (the verdict), and the file records "verified_by".
The own-system point check stays in the raw record ("nums") as a secondary note.
"""
import argparse
import collections
import hashlib
import json
import os
import re
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import adapters as A  # noqa: E402
import crosscheck as X  # noqa: E402
import shootout as SH  # noqa: E402
import latex_report as L  # noqa: E402

SYSTEM = {"FriCAS": "fricas", "Axiom": "axiom", "Maxima": "maxima", "sympy": "sympy", "Wolfram": "wolfram"}


def system_of(lab):
    return "mathilda" if lab.startswith("Mathilda:") else SYSTEM[lab]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--chunk", type=int, default=10)
    a = ap.parse_args()
    bench = os.path.abspath(a.bench)
    rdir = os.path.join(bench, "results")
    date = sorted(f[:-5] for f in os.listdir(rdir) if re.match(r"\d{4}-\d\d-\d\d\.json$", f))[-1]
    rpath = os.path.join(rdir, date + ".json")
    saved = json.load(open(rpath))
    rows = saved["rows"]
    cases = {c["id"]: c for c in SH.load_corpus(os.path.join(bench, "corpus.tsv"))}
    # full answer texts (and the "cut" mark) from the run cache, on copies of the raw records
    full = [{"id": r["id"], "raw": {k: dict(v) for k, v in r["raw"].items()}} for r in rows]
    L.full_texts(bench, full)
    fulltext = {(f["id"], lab): v for f in full for lab, v in f["raw"].items()}

    ver = os.popen(f"{A.MATHILDA} -v").read().split("\n")[0]
    stamp = SH.mathilda_stamp()      # the binary and its .m packages: a rebuild re-judges
    pzq_path = os.path.join(rdir, "pzq.jsonl")
    known = {}
    if os.path.exists(pzq_path):
        for ln in open(pzq_path):
            try:
                rec = json.loads(ln)
                known[rec["key"]] = rec["verdict"]
            except (ValueError, KeyError):
                pass

    verdict, todo = {}, []
    for r in rows:
        c = cases.get(r["id"])
        for lab, raw in r["raw"].items():
            if raw.get("status") != "ANSWER" or c is None:
                continue
            fr = fulltext[(r["id"], lab)]
            text = fr.get("text") or ""
            if fr.get("cut"):
                verdict[(r["id"], lab)] = "cut"
                continue
            try:
                wl = text if system_of(lab) in ("mathilda", "wolfram") else X.to_wl(system_of(lab), text)
            except Exception:  # noqa: BLE001 -- any translator failure: not judged, so unverified
                verdict[(r["id"], lab)] = "untranslatable"
                continue
            if "\n" in wl or '"' in wl:
                verdict[(r["id"], lab)] = "untranslatable"
                continue
            key = hashlib.sha1(f"{stamp}|{c['f']}|{wl}".encode()).hexdigest()
            if key in known:
                verdict[(r["id"], lab)] = known[key]
            else:
                todo.append(((r["id"], lab), key, c["f"], wl))
    print(f"{os.path.basename(bench)}: {len(verdict) + len(todo)} answers, {len(todo)} to judge "
          f"by PossibleZeroQ ({len(verdict)} known or not judgeable)", flush=True)

    tmpdir = tempfile.mkdtemp(prefix="pzq-")
    chunks = [todo[i:i + a.chunk] for i in range(0, len(todo), a.chunk)]
    lock_out = open(pzq_path, "a")

    def work(chunk):
        got = X.run_batch([(f"p{j}", t[2], t[3]) for j, t in enumerate(chunk)], tmpdir)
        return [(t, got.get(f"p{j}", "undecided")) for j, t in enumerate(chunk)]

    done = 0
    with ThreadPoolExecutor(a.jobs) as ex:
        for res in ex.map(work, chunks):
            for (rk, key, _f, _wl), v in res:
                verdict[rk] = v
                lock_out.write(json.dumps({"key": key, "verdict": v, "id": rk[0], "run": rk[1]}) + "\n")
            lock_out.flush()
            done += len(res)
            if done % 500 < a.chunk:
                print(f"  {done}/{len(todo)}", flush=True)
    lock_out.close()

    tally = collections.Counter()
    for r in rows:
        for lab, raw in r["raw"].items():
            v = verdict.get((r["id"], lab))
            if v is None:
                continue
            fr = fulltext[(r["id"], lab)]
            judged = dict(raw, text=fr.get("text") or raw.get("text") or "")
            # the PossibleZeroQ verdict decides: zero -> verified, nonzero -> wrong, else unverified
            r["runs"][lab] = SH.classify_run(system_of(lab), judged, v if v in ("zero", "nonzero") else "undecided")
            raw["pzq"] = v
            tally[(lab, r["runs"][lab])] += 1
    saved["verified_by"] = ("PossibleZeroQ[D[F, x] - f] for every answer (pzq_verify.py, " + ver + ")")
    json.dump(saved, open(rpath, "w"), indent=0)
    by = collections.defaultdict(dict)
    for (lab, cls), n in sorted(tally.items()):
        by[lab][cls] = n
    for lab, d in by.items():
        print(f"  {lab}: {d}")


if __name__ == "__main__":
    main()
