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
RATIO = 2
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

# Benchmark 261 (--plan 261): the follow-on corpus, F at least FOUR times f
# (strict: 4 LeafCount[f] < LeafCount[compactForm[F]]), weighted to the mixed
# classes M and B (a radical over a Log/Exp tower, and below one), with fresh
# seeds (101+) and the integrands of the earlier compact corpora excluded.
# Per-run quotas carry ~20% surplus for the content reclassification and the
# duplicate cut.
PLANS = {
    "257": {"bench": "257-compact-integrals", "ratio": 2, "quotas": QUOTAS, "seeds": SEEDS, "attempts": 40},
    "261": {"bench": "261-compact-integrals-x4", "ratio": 4, "attempts": 150,
            "quotas": {"R": 50, "T": 250, "A": 250, "M": 550, "B": 450},
            "seeds": {"R": [(101, "R", 60)],
                      # seeds past the first build (T 107+, A 107+, M 113+, B2 115+) were added
                      # with --extend: dedup and the exclusions left T, A, M and B short of quota
                      "T": [(s, "T", 50) for s in range(101, 115)],
                      "A": [(s, "A", 50) for s in range(101, 110)],
                      "M": [(s, "M", 55) for s in range(101, 115)],
                      "B": [(s, "B", 40) for s in range(101, 105)] + [(s, "B2", 40) for s in range(105, 118)]},
            "cand_systems": ["fricas", "axiom", "maxima"],
            "exclude": ["257-compact-integrals/gen_pool.tsv", "256-integration-shootout/corpus.tsv",
                        "260-radical-bottom/corpus.tsv", "259-transcendental-decision/corpus.tsv"]},
    # Benchmark 262 (--plan 262): 1000 integrals outside the easy classes of
    # easy.m -- not rational, not in Q(x, ((a x + b)/(c x + d))^(1/m)), not
    # Euler-rationalisable (nor rationalisable by a chain of those substitutions),
    # and not closed by Method -> "DerivativeDivides".  Grammars T3/A3/M3/B3
    # (gen_compact.m), F at least twice f (ratio 2 as in 257, fair mode on).
    # M3 finds almost nothing (a radical of a Log/Exp polynomial whose cofactor
    # is not the kernel's derivative is rarely elementary): its seeds run, and
    # whatever they yield is kept.
    # M runs at ratio 1 (F larger than f): an elementary, non-derivative-divides
    # radical over an Exp tower is built on a unit, Log[P + Sqrt[P^2 - c E^h]],
    # whose answer is never twice its integrand (grammar M4 in gen_compact.m)
    "262": {"bench": "262-hard-integrals", "ratio": 2, "attempts": 20, "exclude_easy": True,
            "class_ratio": {"M": 1},
            "quotas": {"R": 0, "T": 330, "A": 337, "M": 215, "B": 118},
            "seeds": {"R": [],
                      "T": [(s, "T3", 30) for s in range(201, 215)],
                      "A": [(s, "A3", 30) for s in range(201, 216)],
                      "M": [(s, "M3", 15) for s in range(201, 205)] + [(s, "M4", 15) for s in range(205, 223)],
                      "B": [(s, "B3", 30) for s in range(201, 215)]},
            "cand_systems": ["fricas", "axiom", "maxima"],
            "exclude": ["261-compact-integrals-x4/gen_pool.tsv", "257-compact-integrals/gen_pool.tsv",
                        "256-integration-shootout/corpus.tsv", "260-radical-bottom/corpus.tsv",
                        "259-transcendental-decision/corpus.tsv"]},
    # Benchmark 263 (--plan 263): more integrals outside the easy classes, to
    # replace the easy cases easy_classify.py found in benchmarks 256-261 (most of
    # 257 and 261).  The 262 grammars on fresh seeds, plus B4 (gen_compact.m): B3
    # repeated itself (110 distinct integrands in 118 cases).  B4 builds the
    # derivative of a product holding a radical of a cubic or quartic under
    # Exp/Log, so its F is smaller than f: B runs at ratio 0 (any elementary F),
    # M at ratio 1 as in 262, T and A at ratio 2.  Every earlier corpus and pool,
    # 262's included, is excluded.
    "263": {"bench": "263-hard-integrals-ii", "ratio": 2, "attempts": 20, "exclude_easy": True,
            "class_ratio": {"M": 1, "B": 0},
            "quotas": {"R": 0, "T": 500, "A": 550, "M": 400, "B": 400},
            "seeds": {"R": [],
                      "T": [(s, "T3", 30) for s in range(301, 331)],
                      "A": [(s, "A3", 30) for s in range(301, 336)],
                      "M": [(s, "M4", 15) for s in range(301, 341)],
                      "B": [(s, "B4", 15) for s in range(301, 336)] + [(s, "B3", 30) for s in range(336, 341)]},
            "cand_systems": ["fricas", "axiom", "maxima"],
            "exclude": ["262-hard-integrals/gen_pool.tsv", "262-hard-integrals/corpus.tsv",
                        "261-compact-integrals-x4/gen_pool.tsv", "257-compact-integrals/gen_pool.tsv",
                        "256-integration-shootout/corpus.tsv", "260-radical-bottom/corpus.tsv",
                        "259-transcendental-decision/corpus.tsv"]},
}


def content_class(cls, f):
    """The class a case belongs to by what its integrand contains: a radical
    with no Log/Exp kernel is algebraic (A), a kernel with no radical is
    transcendental (T), neither is rational (R); M and B keep their class."""
    rad = bool(re.search(r"\^\(-?\d+/\d+\)|Sqrt\[", f))
    ker = bool(re.search(r"Log\[|E\^|Exp\[", f))
    if rad and ker:
        return cls if cls in ("M", "B") else "B"
    return "A" if rad else ("T" if ker else "R")


EXCLUDE_EASY = False


CLASS_RATIO = {}


def ratio_of(cls):
    """the leaf-count ratio for a class: the plan's, unless class_ratio overrides it"""
    return CLASS_RATIO.get(cls, RATIO)


def ratio_ok(cls, lF, lf):
    r = ratio_of(cls)
    return lF >= 2 * lf if r == 2 else lF > r * lf


def run_gen(cls, seed, grammar, quota, timeout=6 * 3600):
    env = dict(os.environ, SEED=str(seed), CLASS=cls, GRAMMAR=grammar, QUOTA=str(quota),
               RATIO=str(ratio_of(cls)), ATTEMPTS=str(ATTEMPTS), EXCLUDE_EASY="1" if EXCLUDE_EASY else "0",
               FAIR="1" if EXCLUDE_EASY else "0")
    p = subprocess.run([A.MATHILDA, "-file", os.path.join(HERE, "gen_compact.m")], env=env,
                       capture_output=True, text=True, timeout=timeout, cwd=A.REPO)
    stat = [ln for ln in p.stderr.splitlines() if ln.startswith("class ")]
    return cls, seed, p.stdout.splitlines(), stat[-1] if stat else p.stderr[-300:]


def attempt_no(case_id):
    m = re.match(r"c[A-Z]-(\d+)-(\d+)$", case_id)
    return (int(m.group(1)), int(m.group(2))) if m else (99, 10 ** 9)


def leaves(text):
    return len(re.findall(r"[A-Za-z_%][A-Za-z0-9_]*|\d+", text))


# The systems that may close a candidate Mathilda did not, in order.  Plan 257
# asked FriCAS alone; plan 261 also asks Axiom and Maxima, so an integrand only
# one of the other systems closes is in the corpus too (family c<class>-<sys>).
CAND_SYSTEMS = ["fricas"]
CAND_RUNNERS = {"fricas": lambda c, t, d: A.run_fricas_one(c, t, d),
                "axiom": lambda c, t, d: A.run_axiom_one(c, t, d),
                "maxima": lambda c, t, d: A.run_maxima_one(c, t, d)}


def fricas_check(line, tmpdir):
    """The first system of CAND_SYSTEMS that closes the candidate with a
    verified, elementary, large enough answer makes it a case."""
    p = line.split("\t")
    pts = p[7].strip("{}")
    c = {"id": p[1], "class": p[2], "fullform": p[9], "tree": T.parse(p[9]),
         "points_list": [s.strip() for s in pts.split(",")]}
    lf = int(p[10]) if len(p) > 10 else 10 ** 6
    for system in CAND_SYSTEMS:
        if not A.available(system):
            continue
        try:
            rec = CAND_RUNNERS[system](c, 30, tmpdir)
        except Exception:
            continue
        ok = (rec.get("status") == "ANSWER" and rec.get("nums") and all(n == "ok" for n in rec["nums"])
              and not T.special_q(system, rec.get("text", "")))
        nl = leaves(rec.get("text", ""))
        if ok and ratio_ok(p[2][0] if p[2] else "", nl, lf):
            q = p[:10]
            q[0] = "CASE"
            q[3] = "c" + p[2] + "-" + system
            return "\t".join(q)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--plan", default="257", choices=sorted(PLANS),
                    help="which corpus: 257 (ratio 2), 261 (ratio 4, mixed-weighted), 262 or 263 (outside the easy classes)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--timeout", type=float, default=6.0, help="hours per generator run")
    ap.add_argument("--extend", action="store_true",
                    help="reuse the cases of the existing corpus and gen_stats.json; run only the seeds not yet run")
    ap.add_argument("--rerun", default="",
                    help="with --extend: CLS:SEED,... runs to redo (their old cases are dropped first)")
    a = ap.parse_args()
    global BENCH, POOL, QUOTAS, SEEDS, RATIO, ATTEMPTS, EXCLUDE_EASY, CLASS_RATIO
    plan = PLANS[a.plan]
    EXCLUDE_EASY = plan.get("exclude_easy", False)
    BENCH = os.path.join(A.REPO, "benchmarks", plan["bench"])
    POOL = os.path.join(BENCH, "gen_pool.tsv")
    QUOTAS, SEEDS, RATIO, ATTEMPTS = plan["quotas"], plan["seeds"], plan["ratio"], plan["attempts"]
    CLASS_RATIO = plan.get("class_ratio", {})
    CAND_SYSTEMS[:] = plan.get("cand_systems", ["fricas"])
    if a.out is None:
        a.out = os.path.join(BENCH, "corpus.tsv")
    excluded = set()
    for rel in plan.get("exclude", []):
        path = os.path.join(A.REPO, "benchmarks", rel)
        if os.path.exists(path):
            for ln in open(path).read().splitlines():
                q = ln.split("\t")
                if len(q) > 5 and q[0] == "CASE":
                    excluded.add(q[5])
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
    # every finished generator run is kept in parts/<class>-<seed>.tsv (its
    # summary on the first line): an interrupted build resumes, and its cases
    # can be read while the others still run
    parts = os.path.join(BENCH, "parts")
    os.makedirs(parts, exist_ok=True)

    def take(c, s, lines, stat):
        stats["generators"].append({"class": c, "seed": s, "summary": stat})
        print(f"  {c} seed {s}: {stat}", flush=True)
        for ln in lines:
            if ln.startswith("CASE\t"):
                cases[c].append(ln)
            elif ln.startswith("CAND\t"):
                candl[c].append(ln)

    todo = []
    for c in QUOTAS:
        for (sd, g, q) in SEEDS[c]:
            if (c, sd) in done:
                continue
            pf = os.path.join(parts, f"{c}-{sd}.tsv")
            if os.path.exists(pf) and not a.extend:
                lines = open(pf).read().splitlines()
                take(c, sd, lines[1:], lines[0] + " (resumed)")
            else:
                todo.append((c, sd, g, q))
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(run_gen, c, sd, g, q, int(a.timeout * 3600)) for (c, sd, g, q) in todo]
        for fu in cf.as_completed(futs):
            c, sd, lines, stat = fu.result()
            with open(os.path.join(parts, f"{c}-{sd}.tsv"), "w") as fh:
                fh.write(stat + "\n" + "\n".join(lines) + "\n")
            take(c, sd, lines, stat)
    # --extend keeps the candidates of the earlier runs (already checked, not re-run)
    cpath = os.path.join(BENCH, "candidates.tsv")
    oldc = open(cpath).read().splitlines() if a.extend and os.path.exists(cpath) else []
    with open(cpath, "w") as fh:
        fh.write("\n".join([ln for ln in oldc if ln] + [ln for c in QUOTAS for ln in candl[c]]) + "\n")
    with tempfile.TemporaryDirectory() as tmpdir, cf.ThreadPoolExecutor(a.jobs) as ex:
        for c in QUOTAS:
            got = [r for r in ex.map(lambda ln: fricas_check(ln, tmpdir), candl[c]) if r]
            prev = stats["fricas"].get(c, {"candidates": 0, "accepted": 0})
            if a.extend:
                prev = json.load(open(os.path.join(BENCH, "gen_stats.json")))["fricas"].get(c, prev)
            stats["fricas"][c] = {"candidates": prev["candidates"] + len(candl[c]), "accepted": prev["accepted"] + len(got)}
            print(f"  {c}: {'/'.join(CAND_SYSTEMS)} accepted {len(got)} of {len(candl[c])} candidates", flush=True)
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
        keep = [ln for ln in keep if ln.split("\t")[5] not in excluded]
        keep = [ln for ln in keep if not (ln.split("\t")[5] in seen or seen.add(ln.split("\t")[5]))]
        stats.setdefault("kept", {})[c] = {"available": len(keep), "kept": min(len(keep), QUOTAS[c]),
                                            "fricas_only": sum(1 for ln in keep[:QUOTAS[c]] if "-fricas\t" in ln)}
        out.extend(keep[:QUOTAS[c]])
        pool.extend(keep)
    with open(a.out, "w") as fh:
        fh.write("\n".join(out) + "\n")
    with open(POOL, "w") as fh:
        fh.write("\n".join(pool) + "\n")
    stats["plan"] = {"name": a.plan, "ratio": RATIO, "attempts_per_case": ATTEMPTS, "excluded_integrands": len(excluded)}
    json.dump(stats, open(os.path.join(BENCH, "gen_stats.json"), "w"), indent=1)
    print(f"wrote {len(out)} cases to {a.out}")


if __name__ == "__main__":
    main()
