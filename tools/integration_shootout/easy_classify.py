#!/usr/bin/env python3
"""easy_classify.py -- record which integrals of a benchmark are "easy", and read that record.

    python3 tools/integration_shootout/easy_classify.py --bench benchmarks/261-compact-integrals-x4 [--jobs 8]

Runs easy.m's easyClass on every integrand of <bench>/corpus.tsv and writes
<bench>/easy.tsv: one line "id <TAB> class" per case, the class being one of
rational, degenerate, logradical, linrad(-sub), euler(-sub), nested(-sub), dd or
none (see easy.m).  Every case that is not "none" is an easy integral: a rational
function, one a substitution rationalises, or one DerivativeDivides closes.  The
shootout, the FriCAS comparison and the LaTeX report leave those out unless they
are run with --include-easy, so a benchmark measures only integrals that need a
decision procedure.

The classification uses the Mathilda binary of adapters.MATHILDA (the dd class
depends on the DerivativeDivides method of that build); the version is recorded
on the header line of easy.tsv.

easy_ids(bench) is the reader the other tools call.
"""
import argparse
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import adapters as A  # noqa: E402

DRIVER = r"""
Get["@EASYM@"];
lines = Select[StringSplit[ReadString["@IN@"], "\n"], StringLength[#] > 0 &];
Scan[Function[ln, Module[{p = StringSplit[ln, "\t"], c},
  c = CheckAbort[easyClass[ToExpression[p[[2]]]], "abort"];
  Print["@@E\t", p[[1]], "\t", c]]], lines];
"""


def easy_ids(bench):
    """ids of the cases of `bench` that easy.tsv classes as easy (empty when there is no easy.tsv)"""
    path = os.path.join(bench, "easy.tsv")
    if not os.path.exists(path):
        return set()
    out = set()
    for ln in open(path):
        if ln.startswith("#"):
            continue
        p = ln.rstrip("\n").split("\t")
        if len(p) == 2 and p[1] not in ("none", "abort"):   # abort: undecided, kept
            out.add(p[0])
    return out


def has_easy_record(bench):
    return os.path.exists(os.path.join(bench, "easy.tsv"))


def classify(items, jobs, chunk=40):
    """[(id, f)] -> {id: class}; a crashed chunk is rerun one case at a time"""
    tmp = tempfile.mkdtemp(prefix="easy-")

    def run(part):
        fin = os.path.join(tmp, f"in{id(part)}.tsv")
        fm = os.path.join(tmp, f"d{id(part)}.m")
        with open(fin, "w") as fh:
            fh.writelines(f"{i}\t{f}\n" for i, f in part)
        with open(fm, "w") as fh:
            fh.write(DRIVER.replace("@EASYM@", os.path.join(HERE, "easy.m")).replace("@IN@", fin))
        # easyClass is at most ~70 s per case (dd: 30 s Integrate + 30 s PossibleZeroQ)
        out, _err, _to, _w = A.run_proc([A.MATHILDA, "-file", fm], 60 + 90 * len(part))
        got = {}
        for ln in out.splitlines():
            p = ln.split("\t")
            if len(p) == 3 and p[0] == "@@E":
                got[p[1]] = p[2]
        return got

    parts = [items[i:i + chunk] for i in range(0, len(items), chunk)]
    res = {}
    with ThreadPoolExecutor(jobs) as ex:
        for got in ex.map(run, parts):
            res.update(got)
    missing = [it for it in items if it[0] not in res]
    if missing:
        with ThreadPoolExecutor(jobs) as ex:
            for it, got in zip(missing, ex.map(lambda it: run([it]), missing)):
                res[it[0]] = got.get(it[0], "abort")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    bench = os.path.abspath(a.bench)
    items = []
    for ln in open(os.path.join(bench, "corpus.tsv")):
        p = ln.rstrip("\n").split("\t")
        if len(p) > 5 and p[0] == "CASE":
            items.append((p[1], p[5]))
    res = classify(items, a.jobs)
    ver = subprocess.run([A.MATHILDA, "-v"], capture_output=True, text=True).stdout.splitlines()[0]
    with open(os.path.join(bench, "easy.tsv"), "w") as fh:
        fh.write(f"# easy.m classes of corpus.tsv ({ver}); every class but none is left out of the shootout\n")
        fh.writelines(f"{i}\t{res[i]}\n" for i, _ in items)
    n_easy = sum(1 for i, _ in items if res[i] not in ("none", "abort"))
    print(f"{os.path.basename(bench)}: {len(items)} cases, {n_easy} easy, {len(items) - n_easy} kept")


if __name__ == "__main__":
    main()
