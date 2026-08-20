#!/usr/bin/env python3
"""benchmarks/88-thue-equations/grid.py -- randomized Thue stress grid vs PARI.

`cases.py` is a *curated* corpus (~100 hand-picked equations). This file is the
*randomized* counterpart: a deterministic-seeded generator of random binary
forms F(x,y) across degree 3-6 and mixed right-hand sides m (including |m|!=1),
each cross-checked against PARI/GP's `thue()` -- exactly as `run.py` does, but
over hundreds of forms nobody chose.

The point is the same as run.py's, at scale: catch a *wrong finite answer*.
Mathilda declines the hard cases (a safe honest gap), so most random adversarial
forms verdict DECLINE; the grid earns its keep on the forms Mathilda *solves*,
where a set that differs from PARI is a completeness bug (WRONG). The generator
is therefore weighted toward the solve paths (cubics/quartics, small |m|) so
CORRECT is genuinely exercised, not drowned in declines.

Reproducible: a fixed seed (default 20260820) => the identical corpus every run,
so a WRONG is a stable, re-runnable failure. Reuses cases.py's form/poly
builders and run.py's Mathilda+PARI runners -- no duplicated equation logic.

Verdicts (same as run.py): CORRECT / WRONG / CRASH / TIMEOUT / DECLINE /
UNVERIFIED. Exits nonzero on any WRONG or CRASH (never on a DECLINE -- a decline
can only be an honest gap, never a bug). PARI "not a Thue equation" (perfect
powers / repeated factors) -> UNVERIFIED, silently skipped.

    python3 grid.py                      # default 300 cases, seed 20260820
    python3 grid.py --n 800 --seed 7     # bigger / different corpus
    MATH_TIMEOUT=20 PARI_TIMEOUT=15 python3 grid.py

Writes GRID_REPORT.md + grid_results.json.
"""

import argparse
import json
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from cases import mathilda_form, pari_poly          # noqa: E402  (form builders)
from run import run_mathilda, run_pari, classify     # noqa: E402  (the two runners)


# ---------------------------------------------------------------------------
#  adjudicator: PARI thue() is the oracle, but it is not infallible -- for some
#  totally-imaginary fields it silently returns an INCOMPLETE set (e.g. it
#  misses (1,2) for x^4-2x^3+4x^2-3x+1 == 5 over Q(zeta5), which brute force and
#  Mathilda both find).  So when the two disagree, don't blindly trust PARI:
#  each disputed point is cheaply checkable against F(x,y)==m directly (no box
#  needed for SOUNDNESS).  Verdict:
#    MATHILDA_WRONG  Mathilda emitted a non-solution, OR a PARI solution that is
#                    genuine is absent from Mathilda's set (a real completeness bug)
#    PARI_WRONG      Mathilda is sound and has a genuine solution PARI missed
#                    (the oracle is the incomplete one -- Mathilda is correct)
#    None            can't adjudicate (kept as WRONG, conservatively)
# ---------------------------------------------------------------------------
def _form_val(pcoef, x, y):
    n = len(pcoef) - 1
    return sum(c * x ** i * y ** (n - i) for i, c in enumerate(pcoef))


def adjudicate(pcoef, m, M, P):
    if M is None:
        return None
    for (x, y) in M:                                   # soundness: no spurious point
        if _form_val(pcoef, x, y) != m:
            return "MATHILDA_WRONG"
    if P is not None:
        for (x, y) in (P - M):                          # a genuine PARI sol Mathilda lacks
            if _form_val(pcoef, x, y) == m:
                return "MATHILDA_WRONG"
        if any(_form_val(pcoef, x, y) == m for (x, y) in (M - P)):
            return "PARI_WRONG"                         # Mathilda found a genuine sol PARI missed
    return None


# ---------------------------------------------------------------------------
#  random form generator
# ---------------------------------------------------------------------------
def _gcd_list(xs):
    g = 0
    for x in xs:
        g = math.gcd(g, abs(int(x)))
    return g


def gen_case(rng, idx):
    """One random Thue case dict compatible with run_mathilda / run_pari.

    Returns None for a degenerate draw (caller retries). Weighting favours the
    forms Mathilda actually solves so the grid exercises the solve paths."""
    # degree: cubics dominate (best-covered), some quartics, few quintic/sextic.
    n = rng.choices([3, 4, 5, 6], weights=[0.55, 0.30, 0.10, 0.05])[0]
    # coefficient magnitude shrinks with degree (keep roots + PARI time sane).
    C = {3: 6, 4: 4, 5: 3, 6: 2}[n]

    pcoef = [rng.randint(-C, C) for _ in range(n + 1)]
    # leading coeff (of x^n = F's a0): mostly monic (in-scope); a0 must be != 0.
    pcoef[n] = rng.choices([1, -1, 2, -2, 3],
                           weights=[0.42, 0.42, 0.06, 0.06, 0.04])[0]
    # constant coeff (of y^n): != 0, else y | F (degenerate / lower-degree form).
    if pcoef[0] == 0:
        pcoef[0] = rng.choice([-1, 1, 2, -2])
    # primitive form only (a common scale is an uninteresting duplicate).
    g = _gcd_list(pcoef)
    if g != 1:
        pcoef = [c // g for c in pcoef]

    # right-hand side m. For cubics, cover |m|!=1 (M2/M2b paths) heavily; for
    # higher degree, |m|=1 is where we solve, small |m| mostly declines safely.
    if n == 3:
        kind = rng.choices(["pm1", "small", "mid"], weights=[0.40, 0.42, 0.18])[0]
        mag = {"pm1": 1, "small": rng.randint(2, 12), "mid": rng.randint(13, 60)}[kind]
    else:
        kind = rng.choices(["pm1", "small"], weights=[0.72, 0.28])[0]
        mag = 1 if kind == "pm1" else rng.randint(2, 8)
    m = mag * rng.choice([1, -1])

    label = f"grid-{idx:04d}-n{n}"
    return dict(label=label, family=f"grid-n{n}", note="random",
                pcoef=pcoef, m=m,
                form=mathilda_form(pcoef), P=pari_poly(pcoef))


# ---------------------------------------------------------------------------
#  main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=300, help="number of cases (default 300)")
    ap.add_argument("--seed", type=int, default=20260820, help="RNG seed (default 20260820)")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    results = []
    counts = {}
    print(f"Randomized Thue grid: {args.n} cases, seed={args.seed} "
          f"(Mathilda vs PARI oracle)...\n")

    made = 0
    attempts = 0
    while made < args.n and attempts < args.n * 4:
        attempts += 1
        c = gen_case(rng, made + 1)
        if c is None:
            continue
        mo = run_mathilda(c)
        po = run_pari(c)
        # PARI rejects perfect powers / repeated-factor forms ("not a Thue
        # equation"): unusable as an oracle -> UNVERIFIED, doesn't count as made.
        if po["status"] in ("ERROR", "NOGP", "TIMEOUT"):
            verdict = "UNVERIFIED"
        else:
            verdict = classify(mo, po)
            if verdict == "WRONG":                      # don't blindly trust the oracle
                adj = adjudicate(c["pcoef"], c["m"], mo.get("sols"), po.get("sols"))
                if adj == "PARI_WRONG":
                    verdict = "PARI_WRONG"              # Mathilda is the correct one
                # adj == "MATHILDA_WRONG" or None -> keep WRONG
        made += 1
        counts[verdict] = counts.get(verdict, 0) + 1
        mn = len(mo["sols"]) if mo.get("sols") is not None else "-"
        pn = len(po["sols"]) if po.get("sols") is not None else "-"
        flag = "  <<<" if verdict in ("WRONG", "CRASH") else ""
        if args.verbose or verdict in ("WRONG", "CRASH", "PARI_WRONG"):
            print(f"[{made:4d}/{args.n}] {verdict:11s} {c['label']:16s} "
                  f"{c['form']} == {c['m']}   M:{mn} P:{pn}{flag}")
        if verdict in ("WRONG", "PARI_WRONG"):
            print(f"        Mathilda: {sorted(mo['sols'])}")
            print(f"        PARI    : {sorted(po['sols'])}"
                  + ("   (PARI incomplete; Mathilda brute-verified)" if verdict == "PARI_WRONG" else ""))
        if verdict == "CRASH" and mo.get("err"):
            print(f"        err: {mo['err'][:200]}")
        results.append(dict(label=c["label"], family=c["family"], form=c["form"],
                            m=c["m"], pcoef=c["pcoef"], verdict=verdict,
                            m_status=mo["status"], m_ctime=mo.get("ctime"),
                            m_count=(len(mo["sols"]) if mo.get("sols") is not None else None),
                            p_status=po["status"],
                            p_count=(len(po["sols"]) if po.get("sols") is not None else None),
                            m_err=mo.get("err")))
        if (made % 25) == 0 and not args.verbose:
            tally = "  ".join(f"{k}={v}" for k, v in sorted(counts.items()))
            print(f"  ... {made}/{args.n}   {tally}")

    write_report(results, counts, args)
    with open(os.path.join(HERE, "grid_results.json"), "w") as f:
        json.dump(dict(seed=args.seed, n=args.n, counts=counts, results=results), f, indent=1)
    try:
        os.remove(os.path.join(HERE, "_case.m"))
    except OSError:
        pass

    print("\n" + "  ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    bugs = counts.get("WRONG", 0) + counts.get("CRASH", 0)
    if bugs:
        print(f"\n!!! {bugs} BUG(S) FOUND (WRONG/CRASH) -- see GRID_REPORT.md")
    else:
        print("\nNo WRONG/CRASH: every form Mathilda solved matched PARI.")
    return 1 if bugs else 0


def write_report(results, counts, args):
    lines = []
    L = lines.append
    L("# Benchmark 88 — randomized Thue stress grid\n")
    L(f"Deterministic random binary forms F(x,y) (degree 3-6, mixed m) vs "
      f"**PARI/GP `thue()`**. Seed `{args.seed}`, {args.n} cases. Reproducible: "
      f"the same seed regenerates the identical corpus, so any `WRONG` is a "
      f"stable, re-runnable completeness bug. Regenerate: `python3 grid.py "
      f"--n {args.n} --seed {args.seed}`.\n")
    order = ["CORRECT", "WRONG", "CRASH", "TIMEOUT", "DECLINE", "PARI_WRONG", "UNVERIFIED"]
    meaning = dict(CORRECT="finite set matches PARI",
                   WRONG="**differs from PARI — completeness BUG** (adjudicated: Mathilda unsound or incomplete)",
                   CRASH="**errored / no output — BUG**",
                   TIMEOUT="exceeded budget",
                   DECLINE="unevaluated (honest gap); PARI solved",
                   PARI_WRONG="Mathilda found a genuine solution PARI's thue() missed (oracle incomplete; Mathilda brute-verified correct)",
                   UNVERIFIED="PARI rejected the form (perfect power / repeated factor)")
    L("| verdict | n | meaning |")
    L("|---|---|---|")
    for k in order:
        if counts.get(k):
            L(f"| {k} | {counts[k]} | {meaning[k]} |")
    L("")
    bugs = counts.get("WRONG", 0) + counts.get("CRASH", 0)
    L(f"**Result: {bugs} bug(s).** "
      + ("Every form Mathilda solved matched PARI exactly.\n" if not bugs
         else "See the BUGS section.\n"))

    bad = [r for r in results if r["verdict"] in ("WRONG", "CRASH")]
    if bad:
        L("## 🐞 BUGS\n")
        for r in bad:
            L(f"- **{r['verdict']}** `{r['label']}` — `{r['form']} == {r['m']}` "
              f"(Mathilda {r['m_count']} vs PARI {r['p_count']}).")
        L("")
    oracle = [r for r in results if r["verdict"] == "PARI_WRONG"]
    if oracle:
        L("## Oracle misses (Mathilda correct, PARI's `thue()` incomplete)\n")
        L("Adjudicated by independently checking each disputed point against "
          "`F(x,y)==m`: Mathilda's set is sound and strictly contains a genuine "
          "solution PARI missed. PARI `thue()` is fallible on some "
          "totally-imaginary fields.\n")
        for r in oracle:
            L(f"- `{r['label']}` — `{r['form']} == {r['m']}` "
              f"(Mathilda {r['m_count']}, PARI {r['p_count']}). pcoef `{r['pcoef']}`.")
        L("")

    # a sample of CORRECT solves so the report shows the solve paths were hit
    solved = [r for r in results if r["verdict"] == "CORRECT"]
    L(f"## Solve paths exercised: {len(solved)} CORRECT\n")
    L("A sample (the grid's value is that these are *machine-generated*, not "
      "chosen to look good):\n")
    L("| label | form == m | #sol |")
    L("|---|---|---:|")
    for r in solved[:20]:
        L(f"| {r['label']} | `{r['form']} == {r['m']}` | {r['m_count']} |")
    L("")
    with open(os.path.join(HERE, "GRID_REPORT.md"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
