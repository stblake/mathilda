#!/usr/bin/env python3
"""verify_compact.py -- re-check benchmark 257's corpus.

    python3 tools/integration_shootout/verify_compact.py [corpus.tsv]

For every case whose reference F is Mathilda's (F != None), Mathilda re-checks
PossibleZeroQ[D[F] - f] and LeafCount[compactForm[F]] >= 2 LeafCount[f]
(gen_compact.m's compactForm), and that F is elementary (elementaryQ, the
allowlist of heads the generator accepts).  FriCAS-reference cases (F = None) were checked
by the harness's numeric derivative test when they were accepted.  Also checks
per-case structure: unique ids and integrands, three points, a known class.
Exit 1 on any failure.
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import adapters as A  # noqa: E402

CORPUS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(A.REPO, "benchmarks", "257-compact-integrals", "corpus.tsv")


def main():
    rows = [ln.split("\t") for ln in open(CORPUS).read().splitlines() if ln.strip()]
    bad = []
    ids = [r[1] for r in rows]
    fs = [r[5] for r in rows]
    if len(set(ids)) != len(ids) or len(set(fs)) != len(fs):
        bad.append("duplicate ids or integrands")
    for r in rows:
        if r[2] not in "RTAMB" or len(r[7].strip("{}").split(",")) != 3:
            bad.append(f"{r[1]}: class or points")
    ref = [r for r in rows if r[6] != "None"]
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as fh:
        fh.write(f'Get["{os.path.join(HERE, "gen_compact.m")}"];\n')
        fh.write("chk[id_, f_, F_] := Module[{z = TrueQ[TimeConstrained[PossibleZeroQ[D[F, x] - f], 60, False]],"
                 " lr = LeafCount[TimeConstrained[compactForm[F], 20, F]] >= 2 LeafCount[f], el = elementaryQ[F]},"
                 " If[!(z && lr && el), Print[\"FAIL\\t\", id, \"\\t\", z, \"\\t\", lr, \"\\t\", el]]];\n")
        for r in ref:
            fh.write(f'chk["{r[1]}", {r[5]}, {r[6]}];\n')
        fh.write('Print["DONE"];\n')
        path = fh.name
    env = dict(os.environ, QUOTA="0")
    p = subprocess.run([A.MATHILDA, "-file", path], capture_output=True, text=True, env=env, cwd=A.REPO)
    os.unlink(path)
    out = p.stdout.splitlines()
    bad += [ln for ln in out if ln.startswith("FAIL")]
    if "DONE" not in out:
        bad.append("Mathilda did not finish: " + p.stderr[-300:])
    by = {}
    for r in rows:
        by[r[2]] = by.get(r[2], 0) + 1
    print(f"{len(rows)} cases {by}; {len(ref)} Mathilda-reference cases re-checked; {len(bad)} failures")
    for b in bad:
        print("  " + b)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
