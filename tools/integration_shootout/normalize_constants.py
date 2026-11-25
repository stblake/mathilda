"""normalize_constants.py -- remove the constant factors from a benchmark's integrands.

    python3 tools/integration_shootout/normalize_constants.py --bench benchmarks/264-mixed-integrals [--dry-run]

Integrate[c f, x] is c Integrate[f, x], so a constant factor of an integrand (the -1/2
of -1/2 (x^2 - 3 x - 2)/Sqrt[g], or the Sqrt[2] inside Sqrt[4 + 2 E^(-2 x) + 8 x + 4 x^2])
is no part of the problem.  stripConstantsK (common.m) writes f = k g with g free of
them; the corpus line then carries g, F/k for the reference, and g's values at the same
points.  A case is rewritten only when k is a constant, f - k g is zero and (with a
reference) PossibleZeroQ[D[F/k, x] - g] is True; any other case is reported and kept.

Two integrands that differ only by a constant become the same integral; the later one
is reported as a duplicate, and --drop-duplicates leaves it out of the corpus.
The case ids are kept, and the shootout cache is keyed by the integrand's text, so a
rerun of the shootout runs every system on exactly the rewritten cases.
"""
import argparse
import collections
import concurrent.futures as cf
import os
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
MATHILDA = os.environ.get("SHOOTOUT_MATHILDA", os.path.join(REPO, "Mathilda"))


def run_chunk(lines, tmpdir):
    fd, inp = tempfile.mkstemp(suffix=".tsv", dir=tmpdir)
    with os.fdopen(fd, "w") as fh:
        fh.write("".join(ln + "\n" for ln in lines))
    out = inp + ".out"
    subprocess.run([MATHILDA, "-file", os.path.join(HERE, "normalize_constants.m")], cwd=REPO,
                   env=dict(os.environ, IN=inp, OUT=out), capture_output=True, timeout=3600)
    got = {}
    if os.path.exists(out):
        for ln in open(out):
            p = ln.rstrip("\n").split("\t")
            got[p[1]] = p
    return got


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", required=True)
    ap.add_argument("--dry-run", action="store_true", help="count, change nothing")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--drop-duplicates", action="store_true",
                    help="leave out a case whose rewritten integrand is an earlier case's")
    a = ap.parse_args()
    path = os.path.join(a.bench, "corpus.tsv")
    lines = [ln.rstrip("\n") for ln in open(path)]
    cases = [ln for ln in lines if ln.startswith("CASE\t")]
    tmpdir = tempfile.mkdtemp(prefix="normc-")
    step = max(1, (len(cases) + a.jobs * 4 - 1) // (a.jobs * 4))
    res = {}
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        for got in ex.map(lambda i: run_chunk(cases[i:i + step], tmpdir), range(0, len(cases), step)):
            res.update(got)
    count = collections.Counter()
    new, seen, dups = {}, {}, []
    for ln in cases:
        cid = ln.split("\t")[1]
        r = res.get(cid)
        if r is None:
            count["no result (kept)"] += 1
            f = ln.split("\t")[5]
        elif r[0] == "SAME":
            count["unchanged"] += 1
            f = ln.split("\t")[5]
        elif r[2] != "ok":
            count[f"not rewritten: {r[2]}"] += 1
            f = ln.split("\t")[5]
        else:
            count["rewritten"] += 1
            new[cid] = "\t".join(r[3:])
            f = r[8]
        if f in seen:
            dups.append((cid, seen[f]))
        seen.setdefault(f, cid)
    print(f"{a.bench}: {len(cases)} cases")
    for k, v in sorted(count.items()):
        print(f"  {k}: {v}")
    print(f"  duplicates after the rewrite: {len(dups)}" + ("".join(f"\n    {c} = {o}" for c, o in dups[:20])))
    if a.dry_run or not new:
        return
    drop = {c for c, _ in dups} if a.drop_duplicates else set()
    with open(path, "w") as fh:
        for ln in lines:
            if ln.startswith("CASE\t"):
                cid = ln.split("\t")[1]
                if cid in drop:
                    continue
                ln = new.get(cid, ln)
            fh.write(ln + "\n")
    if drop:
        print(f"  left out {len(drop)} duplicates")
    print(f"  wrote {path}")


if __name__ == "__main__":
    main()
