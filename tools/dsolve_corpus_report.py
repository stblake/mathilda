#!/usr/bin/env python3
"""
dsolve_corpus_report.py -- bucket the DSolve corpus TSV by Maple classification
and rank the gaps (UNEVAL/FAIL) so the campaign can target the biggest first.

Input: the per-case TSV emitted by tests/test_dsolve_corpus.c on stdout, lines
    <P|F|U|S>\t<label>\t<classif>
Usage: python3 tools/dsolve_corpus_report.py corpus.tsv [--examples N]
"""
import sys, re, argparse, collections

# Ordered (name, predicate) -- first match wins.  predicate(cl) on the lowercased
# classification string.
def has(cl, *subs): return any(s in cl for s in subs)
RULES = [
    ("system",              lambda cl: 'system' in cl),
    ("1st_Abel",            lambda cl: '_abel' in cl),
    ("1st_solvable_for_yx", lambda cl: 'y=_g(x' in cl or 'x=_g(y' in cl or '_dalembert' in cl),
    ("1st_Chini",           lambda cl: '_chini' in cl),
    ("1st_Riccati",         lambda cl: '_riccati' in cl),
    ("1st_Bernoulli",       lambda cl: '_bernoulli' in cl),
    ("1st_with_symmetry",   lambda cl: '_1st_order' in cl and '_with_symmetry' in cl),
    ("1st_other",           lambda cl: '_1st_order' in cl),
    ("ortho_poly",          lambda cl: has(cl, '_gegenbauer', '_jacobi', '_laguerre', '_hermite')),
    ("Emden_Fowler",        lambda cl: has(cl, '_emden', '_fowler')),
    ("elliptic",            lambda cl: has(cl, '_ellipsoidal', '_elliptic')),
    ("Lienard",             lambda cl: '_lienard' in cl),
    ("Liouville_2nd",       lambda cl: '_liouville' in cl),
    ("2nd_reducible_mu",    lambda cl: '_2nd_order' in cl and has(cl, '_reducible', '_mu_', '_missing_x', '_missing_y', '_exact')),
    ("2nd_linear",          lambda cl: '_2nd_order' in cl and has(cl, '_with_linear_symmetries', '_linear', '_nonhomogeneous')),
    ("2nd_other",           lambda cl: '_2nd_order' in cl),
    ("3rd_high_reducible",  lambda cl: has(cl, '_3rd_order', '_high_order') and has(cl, '_reducible', '_mu_', '_missing_x', '_missing_y', '_exact')),
    ("3rd_high_linear",     lambda cl: has(cl, '_3rd_order', '_high_order')),
    ("Bessel_special",      lambda cl: has(cl, '_bessel', '_kummer', '_whittaker', '_titchmarsh')),
    ("NONE",                lambda cl: 'none' in cl),
    ("rational_misc",       lambda cl: True),
]

def bucket(cl):
    cl = cl.lower()
    for name, pred in RULES:
        if pred(cl): return name
    return "rational_misc"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv"); ap.add_argument("--examples", type=int, default=4)
    args = ap.parse_args()
    rows = []
    for line in open(args.tsv):
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2: continue
        code = parts[0].strip()
        label = parts[1] if len(parts) > 1 else ""
        cl = parts[2] if len(parts) > 2 else ""
        cl = cl.strip().strip('"')
        rows.append((code, label, cl))

    tot = collections.Counter(r[0] for r in rows)
    scalar = [r for r in rows if r[0] != 'S']
    npass = sum(1 for r in scalar if r[0] == 'P')
    print("== Overall ==")
    print(f"  records: {len(rows)}   PASS {tot['P']}  FAIL {tot['F']}  UNEVAL {tot['U']}  SKIP(sys) {tot['S']}")
    n_scalar = len(scalar)
    print(f"  scalar: {n_scalar}   solved {npass} ({100.0*npass/max(1,n_scalar):.1f}%)   "
          f"gap {n_scalar-npass}")

    b = collections.defaultdict(lambda: collections.Counter())
    ex = collections.defaultdict(list)
    for code, label, cl in scalar:
        bk = bucket(cl); b[bk][code] += 1; b[bk]['tot'] += 1
        if code in ('U', 'F') and len(ex[bk]) < args.examples:
            ex[bk].append(label)

    order = sorted(b, key=lambda k: -(b[k]['U'] + b[k]['F']))
    print("\n== Scalar buckets, ranked by gap (UNEVAL+FAIL) ==")
    print(f"  {'bucket':22s} {'tot':>4s} {'PASS':>5s} {'UNEV':>5s} {'FAIL':>5s} {'pass%':>6s}   examples")
    for k in order:
        c = b[k]; t = c['tot']; p = c['P']
        gap = c['U'] + c['F']
        star = ' <=GAP' if gap else ''
        print(f"  {k:22s} {t:4d} {p:5d} {c['U']:5d} {c['F']:5d} {100.0*p/max(1,t):6.1f}   "
              f"{', '.join(ex[k][:args.examples])}{star}")

    fails = [r for r in scalar if r[0] == 'F']
    if fails:
        print(f"\n== {len(fails)} FAIL (wrong closed form -- investigate) ==")
        for _, label, cl in fails[:40]:
            print(f"   {label}\t{cl[:70]}")

if __name__ == "__main__":
    main()
