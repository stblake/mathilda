#!/usr/bin/env python3
"""stress.py -- memory-safe stress driver for the Bronstein package.

Generates seeded batches (gen.m), runs each batch in its own Mathilda process
(run.m) under a wall-clock timeout, and aggregates RESULT lines by class,
family and decline reason.  A FALSE_CERT, WRONG or CRASH anywhere exits 1;
--max-open N also fails when more than N cases end DECLINE/TIMEOUT/UNEVAL (a
coverage ratchet).  --family random generates the mostly non-elementary random
family, and --fricas cross-checks every "not elementary" verdict and every
decline against FriCAS (fricas_check.py): a DISAGREE exits 1.  Exits 77 (ctest
SKIP) when ./Mathilda, or FriCAS for --fricas, is not installed.

    python3 tools/bronstein_stress/stress.py --seeds 1-10 --count 20 [--method Automatic]
    python3 tools/bronstein_stress/stress.py --family random --seeds 1-8 --fricas
"""
import argparse, collections, os, subprocess, sys, tempfile, json

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MATH = os.path.join(ROOT, "Mathilda")
HERE = os.path.dirname(os.path.abspath(__file__))

def run(cmd, env, timeout):
    try:
        p = subprocess.run(cmd, env=env, cwd=ROOT, capture_output=True, text=True, timeout=timeout)
        return p.stdout
    except subprocess.TimeoutExpired as e:
        return (e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or "")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="1-5")
    ap.add_argument("--count", type=int, default=20)
    ap.add_argument("--method", default="Bronstein")
    ap.add_argument("--budget", type=int, default=30)
    ap.add_argument("--json", default=None)
    ap.add_argument("--family", default=None)
    ap.add_argument("--max-open", type=int, default=None)
    ap.add_argument("--fricas", action="store_true")
    a = ap.parse_args()
    if not os.access(MATH, os.X_OK):
        print("SKIP: ./Mathilda not built"); sys.exit(77)
    if a.fricas:
        sys.path.insert(0, HERE)
        import fricas_check
        if not os.access(fricas_check.FRICAS, os.X_OK):
            print("SKIP: FriCAS not installed at " + fricas_check.FRICAS); sys.exit(77)
    lo, hi = (int(v) for v in a.seeds.split("-"))
    results = []
    for seed in range(lo, hi + 1):
        env = dict(os.environ, SEED=str(seed), COUNT=str(a.count))
        if a.family: env["FAMILY"] = a.family
        cases = [l for l in run([MATH, "-file", os.path.join(HERE, "gen.m")], env, 600).splitlines() if l.startswith("CASE\t")]
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
            fh.write("\n".join(cases) + "\n"); path = fh.name
        env = dict(os.environ, CASES=path, METHOD=a.method, BUDGET=str(a.budget))
        out = run([MATH, "-file", os.path.join(HERE, "run.m")], env, a.budget * (a.count + 2) + 120)
        got = {}
        for l in out.splitlines():
            if l.startswith("RESULT\t"):
                p = l.split("\t") + [""] * 6
                got[p[1]] = dict(id=p[1], fam=p[2], cls=p[3], secs=float(p[4] or 0), detail=p[5])
        for c in cases:
            parts = c.split("\t")
            cid, fam, f = parts[1], parts[2], parts[3]
            r = got.get(cid, dict(id=cid, fam=fam, cls="CRASH", secs=0.0, detail=""))
            r["f"] = f
            results.append(r)
        os.unlink(path)
        print(f"seed {seed}: " + ", ".join(f"{k}={v}" for k, v in collections.Counter(r['cls'] for r in results[-len(cases):]).items()), flush=True)
    cls = collections.Counter(r["cls"] for r in results)
    print("\nTOTAL", len(results), dict(cls))
    byfam = collections.defaultdict(collections.Counter)
    for r in results: byfam[r["fam"]][r["cls"]] += 1
    for f, c in sorted(byfam.items()): print(f"  {f:8s} {dict(c)}")
    why = collections.Counter(r["detail"] for r in results if r["cls"] == "DECLINE")
    print("\nDECLINE reasons:")
    for w, n in why.most_common(): print(f"  {n:4d}  {w}")
    slow = sorted(results, key=lambda r: -r["secs"])[:8]
    print("\nslowest:"); [print(f"  {r['secs']:7.2f}  {r['cls']:9s} {r['fam']:8s} {r['f'][:110]}") for r in slow]
    # on the random family (mostly NOT elementary) a certificate is a claim for
    # FriCAS to check, not a failure; on the known-elementary families it is wrong
    hard = ("WRONG", "CRASH") if a.family == "random" else ("FALSE_CERT", "WRONG", "CRASH")
    bad = [r for r in results if r["cls"] in hard]
    for r in bad: print(f"\n*** {r['cls']} {r['id']} {r['fam']}: {r['f']}\n    {r['detail'][:300]}")
    if a.json: json.dump(results, open(a.json, "w"), indent=1)
    fail = bool(bad)
    nopen = sum(cls[k] for k in ("DECLINE", "TIMEOUT", "UNEVAL"))
    if a.max_open is not None and nopen > a.max_open:
        print(f"\n*** {nopen} open cases (DECLINE/TIMEOUT/UNEVAL) exceed --max-open {a.max_open}"); fail = True
    if a.fricas:
        n = fricas_check.check(results)
        if n["DISAGREE"]: fail = True
    sys.exit(1 if fail else 0)

if __name__ == "__main__":
    main()
