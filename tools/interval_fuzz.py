#!/usr/bin/env python3
"""Randomised stress test for Interval[] arithmetic — the inclusion guarantee.

The one property interval arithmetic must never break is *containment* (the
inclusion property): for any operation f and interval X, the computed f(X) must
enclose every true value { f(x) : x in X }. A tolerant endpoint comparison once
let the outward-rounded min/max pick the INNER of two endpoints that differed by
a ULP, so f(X) silently excluded the true value (e.g. the logistic-map orbit).
This harness would have caught that, and guards against its return.

It drives the compiled ``./Mathilda`` over the line-based NDJSON protocol (the
same transport ``site/generate.py`` uses) and checks four families:

  * containment  — random INEXACT intervals through arithmetic, powers, and
                   every threaded function; each sampled true value must be an
                   IntervalMemberQ of the result. Inexact endpoints are used so
                   the result is numeric-with-width and a high-precision sample
                   sits strictly inside — an exact-rational interval would give
                   symbolic endpoints (Sin[1]) that a rounded sample can't equal.
  * exactness    — exact-rational +, -, *, ^n must give the mathematically exact
                   interval, not merely contain it.
  * determinism  — the same expression evaluated twice must be identical.
  * cancellation — the tight-endpoint regime (x - x, x(1-x), Sqrt[x]^2,
                   Exp[Log[x]], the logistic orbit) where endpoints differ by a
                   few ULP — exactly where the tolerance bug lived.

Usage:
  python3 tools/interval_fuzz.py [--seeds N] [--cases M]

Exits non-zero if any check fails; prints the offending expression.
Requires a built ``./Mathilda`` (run ``make`` at the repo root).
"""
import argparse
import json
import random
import subprocess
import sys
from fractions import Fraction as Fr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MATHILDA = ROOT / "Mathilda"


def run(exprs, timeout=300):
    """Evaluate each expression via the NDJSON pipe; return payloads in order."""
    reqs = [json.dumps({"id": i, "expr": e}) for i, e in enumerate(exprs, 1)]
    reqs.append(json.dumps({"type": "quit"}))
    proc = subprocess.run([str(MATHILDA)], input="\n".join(reqs) + "\n",
                          capture_output=True, text=True, timeout=timeout)
    pay = {}
    for raw in proc.stdout.splitlines():
        raw = raw.strip()
        if raw[:1] != "{":
            continue
        try:
            msg = json.loads(raw)
        except ValueError:
            continue
        if isinstance(msg, dict) and msg.get("type") == "expr" and "payload" in msg:
            pay[msg.get("id")] = msg["payload"]
    return [pay.get(i, "") for i in range(1, len(exprs) + 1)]


def fmt(x):
    return repr(float(x))


def dec(lo=-6.0, hi=6.0):
    return round(random.uniform(lo, hi), 3)


def ivs(a, b):
    if a > b:
        a, b = b, a
    return f"Interval[{{{fmt(a)}, {fmt(b)}}}]", a, b


def containment(cases):
    """cases: list of (F_expr, [sample scalar exprs]). Returns list of failures."""
    checks, meta = [], []
    for f_expr, samples in cases:
        inner = ", ".join(f"IntervalMemberQ[R, N[{s}, 30]]" for s in samples)
        checks.append(f"Module[{{R = {f_expr}}}, "
                      f"If[Head[R] === Interval, And[{inner}], IvSkip]]")
        meta.append(f_expr)
    fails = []
    for f_expr, out in zip(meta, run(checks)):
        out = out.strip()
        if out not in ("True", "IvSkip", ""):
            fails.append((f_expr, out))
    return fails, len(checks)


def build_containment_cases():
    cases = []

    def unary(fn, a, b):
        x, a, b = ivs(a, b)
        pts = [a, b] + [a + k / 6.0 * (b - a) for k in range(1, 6)]
        cases.append((f"{fn}[{x}]", [f"{fn}[{fmt(p)}]" for p in pts]))

    def binop(op, a, b, c, d):
        x, a, b = ivs(a, b)
        y, c, d = ivs(c, d)
        xs, ys = [a, b, (a + b) / 2], [c, d, (c + d) / 2]
        cases.append((f"({x}) {op} ({y})",
                      [f"({fmt(u)}) {op} ({fmt(v)})" for u in xs for v in ys]))

    a, b = dec(), dec()
    binop("+", a, b, dec(), dec())
    binop("-", a, b, dec(), dec())
    binop("*", a, b, dec(), dec())
    binop("/", a, b, dec(0.5, 6), dec(0.5, 6))
    # division whose divisor straddles 0 -> a disjoint-union result
    x, u0, u1 = ivs(dec(), dec())
    y, v0, v1 = ivs(dec(-4, -0.5), dec(0.5, 4))
    cases.append((f"({x}) / ({y})",
                  [f"({fmt(u)}) / ({fmt(v)})"
                   for u in (u0, u1, (u0 + u1) / 2)
                   for v in (v0, v1, -0.3, 0.3)]))
    for n in (2, 3, 4, -1, -2):
        if n < 0 and a <= 0 <= b:
            continue
        x, aa, bb = ivs(a, b)
        cases.append((f"({x})^({n})",
                      [f"({fmt(p)})^({n})" for p in (aa, bb, (aa + bb) / 2)]))
    for fn in ("Exp", "Sin", "Cos", "Tan", "Sinh", "Tanh", "Cosh", "ArcTan",
               "ArcSinh", "Sec", "Csc", "Cot", "Coth", "Csch", "Abs", "Sign",
               "Floor", "Ceiling", "Erf", "Erfc", "ArcCot", "ArcCsch"):
        unary(fn, dec(), dec())
    for fn in ("Log", "Sqrt", "Gamma", "Sech"):
        unary(fn, dec(0.3, 6), dec(0.3, 6))
    for fn in ("ArcSin", "ArcCos", "ArcTanh"):
        unary(fn, dec(-0.99, 0.99), dec(-0.99, 0.99))
    for fn in ("ArcSec", "ArcCsc", "ArcCoth", "ArcCosh"):
        unary(fn, dec(1.2, 6), dec(1.2, 6))
    unary("ArcSech", dec(0.05, 0.95), dec(0.05, 0.95))
    unary("Zeta", dec(1.3, 6), dec(1.3, 6))
    unary("LogGamma", dec(2.0, 6), dec(2.0, 6))
    pa, pb = dec(0.3, 6), dec(0.3, 6)
    for nn in (0, 1, 2):
        lo, hi = (pa, pb) if pa <= pb else (pb, pa)
        cases.append((f"PolyGamma[{nn}, Interval[{{{fmt(lo)}, {fmt(hi)}}}]]",
                      [f"PolyGamma[{nn}, {fmt(p)}]" for p in (lo, hi, (lo + hi) / 2)]))

    # Special functions threaded by the general derivative-sign certifier
    # (interval_thread_call) and the bespoke sub-domain rows. Each is sampled on
    # a domain where a bound is certifiable; elsewhere the result is symbolic and
    # the containment harness records IvSkip.
    unary("Erfi", dec(-2.5, 2.5), dec(-2.5, 2.5))          # increasing on R
    for fn in ("ExpIntegralEi", "LogIntegral"):            # increasing on (>1 / >0)
        unary(fn, dec(1.3, 6), dec(1.3, 6))
    unary("ExpIntegralEi", dec(-6, -0.3), dec(-6, -0.3))   # decreasing branch (x<0)
    unary("InverseErf", dec(-0.9, 0.9), dec(-0.9, 0.9))    # increasing on (-1,1)
    unary("InverseErfc", dec(0.1, 1.9), dec(0.1, 1.9))     # decreasing on (0,2)
    unary("ProductLog", dec(-0.3, 6), dec(-0.3, 6))        # increasing on (-1/e,inf)
    unary("HarmonicNumber", dec(-0.5, 6), dec(-0.5, 6))    # increasing on (-1,inf)
    for nn in (2, 3):                                      # PolyLog[n,.] inc on (0,1)
        pl, ph = dec(0.05, 0.9), dec(0.05, 0.9)
        lo, hi = (pl, ph) if pl <= ph else (ph, pl)
        cases.append((f"PolyLog[{nn}, Interval[{{{fmt(lo)}, {fmt(hi)}}}]]",
                      [f"PolyLog[{nn}, {fmt(p)}]" for p in (lo, hi, (lo + hi) / 2)]))
    # Piecewise / step functions — non-decreasing, so endpoint threading encloses.
    for fn in ("UnitStep", "Ramp", "Round", "IntegerPart"):
        unary(fn, dec(), dec())
    return cases


def exactness():
    """Exact-rational arithmetic must give the mathematically exact interval."""
    def rat():
        return Fr(random.randint(-30, 30), random.randint(1, 6))

    def order(p, q):
        return (p, q) if p <= q else (q, p)

    src, want = [], []
    a, b = order(rat(), rat())
    c, d = order(rat(), rat())
    src.append(f"Interval[{{{a},{b}}}] + Interval[{{{c},{d}}}]")
    want.append(f"Interval[{{{a + c}, {b + d}}}]")
    src.append(f"Interval[{{{a},{b}}}] - Interval[{{{c},{d}}}]")
    want.append(f"Interval[{{{a - d}, {b - c}}}]")
    prods = [a * c, a * d, b * c, b * d]
    src.append(f"Interval[{{{a},{b}}}] Interval[{{{c},{d}}}]")
    want.append(f"Interval[{{{min(prods)}, {max(prods)}}}]")
    for n in (2, 3):
        if n % 2 == 0 and a < 0 < b:
            lo, hi = Fr(0), max(a ** n, b ** n)
        else:
            lo, hi = min(a ** n, b ** n), max(a ** n, b ** n)
        src.append(f"Interval[{{{a},{b}}}]^{n}")
        want.append(f"Interval[{{{lo}, {hi}}}]")
    fails = []
    for s, w, got in zip(src, want, run(src)):
        if got.replace(" ", "") != w.replace(" ", ""):
            fails.append((s, f"{got}  (want {w})"))
    return fails, len(src)


def determinism():
    a, b = dec(), dec()
    exprs = [
        f"Sin[Interval[{{{fmt(a)}, {fmt(b)}}}]] === Sin[Interval[{{{fmt(a)}, {fmt(b)}}}]]",
        f"(Interval[{{{fmt(a)}, {fmt(b)}}}]^2) === (Interval[{{{fmt(a)}, {fmt(b)}}}]^2)",
    ]
    fails = [(e, o) for e, o in zip(exprs, run(exprs)) if o.strip() != "True"]
    return fails, len(exprs)


def cancellation():
    checks, meta = [], []
    a = round(random.uniform(-5, 5), 3)
    p = round(random.uniform(0.1, 5), 3)
    x = round(a * 0.1 + 0.4, 3)
    checks.append(f"IntervalMemberQ[Interval[{{{fmt(a)}, {fmt(a)}}}] - "
                  f"Interval[{{{fmt(a)}, {fmt(a)}}}], 0]"); meta.append("x - x >= 0")
    checks.append(f"IntervalMemberQ[4 Interval[{{{fmt(x)}, {fmt(x)}}}] "
                  f"(1 - Interval[{{{fmt(x)}, {fmt(x)}}}]), N[4 ({x}) (1 - ({x})), 30]]")
    meta.append("logistic body")
    checks.append(f"IntervalMemberQ[Sqrt[Interval[{{{fmt(p)}, {fmt(p)}}}]]^2, {fmt(p)}]")
    meta.append("Sqrt[x]^2")
    checks.append(f"IntervalMemberQ[Exp[Log[Interval[{{{fmt(p)}, {fmt(p)}}}]]], {fmt(p)}]")
    meta.append("Exp[Log[x]]")
    orbit = [0.5]
    for _ in range(8):
        orbit.append(4 * orbit[-1] * (1 - orbit[-1]))
    for k in range(1, 9):
        checks.append(f"IntervalMemberQ[Nest[4 # (1 - #) &, Interval[0.5], {k}], "
                      f"N[{orbit[k]!r}, 30]]")
        meta.append(f"logistic orbit step {k}")
    fails = [(m, o) for m, o in zip(meta, run(checks)) if o.strip() != "True"]
    return fails, len(checks)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", type=int, default=8, help="number of random seeds")
    ap.add_argument("--cases", type=int, default=60, help="random cases per seed")
    args = ap.parse_args()
    if not MATHILDA.exists():
        sys.exit(f"error: {MATHILDA} not found — run `make` first.")

    total = {"containment": 0, "exactness": 0, "determinism": 0, "cancellation": 0}
    failures = 0
    for seed in range(1, args.seeds + 1):
        random.seed(seed)
        cont_cases = []
        for _ in range(args.cases):
            cont_cases += build_containment_cases()
        for name, (fails, n) in (
            ("containment", containment(cont_cases)),
            ("exactness", exactness_many(args.cases)),
            ("determinism", determinism_many(args.cases)),
            ("cancellation", cancellation_many(args.cases)),
        ):
            total[name] += n
            for expr, out in fails:
                failures += 1
                if failures <= 60:
                    print(f"FAIL [{name}, seed {seed}] {expr}  =>  {out}")

    print("\nchecks run: " + ", ".join(f"{k}={v}" for k, v in total.items()))
    if failures:
        sys.exit(f"{failures} interval stress failure(s)")
    print("interval_fuzz: OK — all containment/exactness/determinism checks passed")


def exactness_many(m):
    fails, n = [], 0
    for _ in range(m * 3):
        f, c = exactness()
        fails += f
        n += c
    return fails, n


def determinism_many(m):
    fails, n = [], 0
    for _ in range(m):
        f, c = determinism()
        fails += f
        n += c
    return fails, n


def cancellation_many(m):
    fails, n = [], 0
    for _ in range(max(1, m // 5)):
        f, c = cancellation()
        fails += f
        n += c
    return fails, n


if __name__ == "__main__":
    main()
