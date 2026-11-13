#!/usr/bin/env python3
"""gen_compact.py -- build benchmark 257 (compact integrands, large integrals).

    python3 tools/integration_shootout/gen_compact.py [--jobs 6] [--out PATH]

Runs tools/integration_shootout/gen_compact.m once per (class, seed) to search
compact integrands (see its header).  Its CASE lines are accepted cases with
Mathilda's verified antiderivative; its CAND lines are candidates Mathilda
neither closed nor certified.  Each CAND goes to FriCAS: when FriCAS returns an
answer that passes the harness's own numeric derivative check (adapters.py, the
same check the shootout applies), is free of special functions, and is at
least twice the integrand's leaf count (the FriCAS answer's leaves estimated by
counting its identifiers and numbers), it becomes a case with F = None and
family c<class>-fricas.  So the corpus is not restricted to integrals Mathilda
can already do.  Wolfram is never consulted: it is a comparison column only.

Per class the first QUOTA accepted cases in (seed, attempt) order are kept.
Writes corpus.tsv and gen_stats.json into the benchmark folder.
"""
import argparse
import concurrent.futures as cf
import json
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import adapters as A  # noqa: E402
import translate as T  # noqa: E402

BENCH = os.path.join(A.REPO, "benchmarks", "257-compact-integrals")
# every accepted case before the per-class quota cut, so --extend sees the
# surplus a from-scratch run would (corpus.tsv is the cut)
POOL = os.path.join(BENCH, "gen_pool.tsv")
QUOTAS = {"R": 100, "T": 150, "A": 300, "M": 250, "B": 200}
# (seed, grammar, per-run quota).  The first runs used each class's own grammar
# with quota/len(seeds); the later B runs use grammar B2 (see gen_compact.m),
# added because B's own grammar drew mostly kernel-free integrands.
SEEDS = {"R": [(1, "R", 100)],
         "T": [(1, "T", 75), (2, "T", 75), (3, "T", 20)],
         "A": [(s, "A", 75) for s in (1, 2, 3, 4)],
         "M": [(s, "M", 63) for s in (1, 2, 3, 4)],
         "B": [(s, "B", 34) for s in (1, 2, 3, 4, 5, 6)] + [(s, "B2", 30) for s in range(7, 13)]}


def content_class(cls, f):
    """The class a case belongs to by what its integrand contains: a radical
    with no Log/Exp kernel is algebraic (A), a kernel with no radical is
    transcendental (T), neither is rational (R); M and B keep their class."""
    rad = bool(re.search(r"\^\(-?\d+/\d+\)|Sqrt\[", f))
    ker = bool(re.search(r"Log\[|E\^|Exp\[", f))
    if rad and ker:
        return cls if cls in ("M", "B") else "B"
    return "A" if rad else ("T" if ker else "R")


def run_gen(cls, seed, grammar, quota, timeout=6 * 3600):
    env = dict(os.environ, SEED=str(seed), CLASS=cls, GRAMMAR=grammar, QUOTA=str(quota))
    p = subprocess.run([A.MATHILDA, "-file", os.path.join(HERE, "gen_compact.m")], env=env,
                       capture_output=True, text=True, timeout=timeout, cwd=A.REPO)
    stat = [ln for ln in p.stderr.splitlines() if ln.startswith("class ")]
    return cls, seed, p.stdout.splitlines(), stat[-1] if stat else p.stderr[-300:]


def attempt_no(case_id):
    m = re.match(r"c[A-Z]-(\d+)-(\d+)$", case_id)
    return (int(m.group(1)), int(m.group(2))) if m else (99, 10 ** 9)


def leaves(text):
    return len(re.findall(r"[A-Za-z_%][A-Za-z0-9_]*|\d+", text))


def fricas_check(line, tmpdir):
    p = line.split("\t")
    pts = p[7].strip("{}")
    c = {"id": p[1], "class": p[2], "fullform": p[9], "tree": T.parse(p[9]),
         "points_list": [s.strip() for s in pts.split(",")]}
    rec = A.run_fricas_one(c, 30, tmpdir)
    ok = (rec.get("status") == "ANSWER" and rec.get("nums") and all(n == "ok" for n in rec["nums"])
          and not T.special_q("fricas", rec.get("text", "")))
    lf = int(p[10]) if len(p) > 10 else 10 ** 6
    if ok and leaves(rec["text"]) >= 2 * lf:
        q = p[:10]
        q[0] = "CASE"
        q[3] = "c" + p[2] + "-fricas"
        return "\t".join(q)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--out", default=os.path.join(BENCH, "corpus.tsv"))
    ap.add_argument("--extend", action="store_true",
                    help="reuse the cases of the existing corpus and gen_stats.json; run only the seeds not yet run")
    ap.add_argument("--rerun", default="",
                    help="with --extend: CLS:SEED,... runs to redo (their old cases are dropped first)")
    a = ap.parse_args()
    os.makedirs(BENCH, exist_ok=True)
    stats = {"generators": [], "fricas": {}}
    cases = {c: [] for c in QUOTAS}
    candl = {c: [] for c in QUOTAS}
    done = set()
    if a.extend:
        old = json.load(open(os.path.join(BENCH, "gen_stats.json")))
        redo = {(c, int(n)) for c, n in (t.split(":") for t in a.rerun.split(",") if t)}
        stats["generators"] = [g for g in old["generators"] if (g["class"], g["seed"]) not in redo]
        stats["reruns"] = old.get("reruns", []) + sorted(f"{c}:{n}" for c, n in redo)
        done = {(g["class"], g["seed"]) for g in stats["generators"]}
        pool = POOL if os.path.exists(POOL) else a.out
        for ln in open(pool).read().splitlines():
            m = re.match(r"c([A-Z])-(\d+)-", ln.split("\t")[1])
            if (m.group(1), int(m.group(2))) not in redo:
                cases[m.group(1)].append(ln)
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(run_gen, c, s, g, q) for c in QUOTAS for (s, g, q) in SEEDS[c] if (c, s) not in done]
        for fu in cf.as_completed(futs):
            c, s, lines, stat = fu.result()
            stats["generators"].append({"class": c, "seed": s, "summary": stat})
            print(f"  {c} seed {s}: {stat}", flush=True)
            for ln in lines:
                if ln.startswith("CASE\t"):
                    cases[c].append(ln)
                elif ln.startswith("CAND\t"):
                    candl[c].append(ln)
    with tempfile.TemporaryDirectory() as tmpdir, cf.ThreadPoolExecutor(a.jobs) as ex:
        for c in QUOTAS:
            got = [r for r in ex.map(lambda ln: fricas_check(ln, tmpdir), candl[c]) if r]
            prev = stats["fricas"].get(c, {"candidates": 0, "accepted": 0})
            if a.extend:
                prev = json.load(open(os.path.join(BENCH, "gen_stats.json")))["fricas"].get(c, prev)
            stats["fricas"][c] = {"candidates": prev["candidates"] + len(candl[c]), "accepted": prev["accepted"] + len(got)}
            print(f"  {c}: FriCAS accepted {len(got)} of {len(candl[c])} candidates", flush=True)
            cases[c].extend(got)
    # reclassify by content: a generated case goes to the class its integrand
    # belongs to (the shootout's method map is per class)
    bycls = {c: [] for c in QUOTAS}
    for c in QUOTAS:
        for ln in cases[c]:
            q = ln.split("\t")
            q[2] = content_class(c, q[5])
            q[3] = "c" + q[2] + q[3][2:]
            bycls[q[2]].append("\t".join(q))
    cases = bycls
    out, seen, pool = [], set(), []
    for c in QUOTAS:
        keep = sorted(cases[c], key=lambda ln: (attempt_no(ln.split("\t")[1]), ln.split("\t")[1]))
        keep = [ln for ln in keep if not (ln.split("\t")[5] in seen or seen.add(ln.split("\t")[5]))]
        stats.setdefault("kept", {})[c] = {"available": len(keep), "kept": min(len(keep), QUOTAS[c]),
                                            "fricas_only": sum(1 for ln in keep[:QUOTAS[c]] if "-fricas\t" in ln)}
        out.extend(keep[:QUOTAS[c]])
        pool.extend(keep)
    with open(a.out, "w") as fh:
        fh.write("\n".join(out) + "\n")
    with open(POOL, "w") as fh:
        fh.write("\n".join(pool) + "\n")
    json.dump(stats, open(os.path.join(BENCH, "gen_stats.json"), "w"), indent=1)
    print(f"wrote {len(out)} cases to {a.out}")


if __name__ == "__main__":
    main()
