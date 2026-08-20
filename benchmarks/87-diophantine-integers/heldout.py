#!/usr/bin/env python3
"""benchmarks/87-diophantine-integers/heldout.py -- the HELD-OUT validation set.

`cases.py` is the *developed-against* benchmark: each row was added in the same
commit as the method that solves it, so a green result there measures "the method
works on the example it was written for", not coverage.  This file is the
opposite discipline: equations drawn from standard references (Mordell, Nagell,
classic problem lists, sympy's own docstrings) that the solver was NOT tuned on,
run COLD, and cross-checked against an independent Python brute-force oracle.

`validate.py` consumes this file.  It never trusts Mathilda's answer: it
brute-forces the same box and flags any disagreement -- above all a `{}` (or a
finite set) that the oracle contradicts, which is a *silent wrong answer* (the
exact class the co-designed benchmark cannot see, and the one that a held-out
run first surfaced for unbounded linear systems).

Each case is a dict:

  label   -- unique id
  math    -- the exact Solve[...] expression Mathilda evaluates
  vars    -- ordered variable names (the tuple order Mathilda emits)
  box     -- {var: (lo, hi)} the ORACLE enumerates (the problem's box for a
             bounded case; a generous validation window for an unbounded one)
  sat     -- lambda d: bool, the equation+constraints as a Python predicate on
             a {var: int} assignment; the oracle's ground truth
  kind    -- 'finite'      : bounded; Mathilda must return exactly the oracle set
             'empty'       : the oracle set is empty; Mathilda must return {}
             'param'       : unbounded family; Mathilda's members (over the box)
                             must equal the oracle set
             'decline_ok'  : unevaluated is acceptable (research-grade); but any
                             concrete answer Mathilda *does* give must be correct
  source  -- provenance / note
"""

from math import gcd, isqrt


def _sq(n):
    return n >= 0 and isqrt(n) ** 2 == n


def _icbrt(n):
    """Integer cube root of n (any sign), or None if n is not a perfect cube."""
    neg = n < 0
    a = abs(n)
    r = round(a ** (1 / 3)) if a else 0
    for c in (r - 1, r, r + 1):
        if c >= 0 and c ** 3 == a:
            return -c if neg else c
    return None


# Equation-aware oracles for cases whose box is too large for a naive grid.
# Each returns a sorted list of solution tuples (in the case's `vars` order)
# lying inside the case's box.

def _oracle_mordell(k, xlo, xhi):
    out = []
    for x in range(xlo, xhi + 1):
        v = x ** 3 + k
        if v >= 0 and _sq(v):
            y = isqrt(v)
            out += ([(x, 0)] if y == 0 else [(x, -y), (x, y)])
    return sorted(out)


def _oracle_thue_x3_m2y3(bound):
    out = []
    for x in range(-bound + 1, bound):
        num = x ** 3 - 1
        if num % 2:
            continue
        y = _icbrt(num // 2)
        if y is not None and abs(y) < bound:
            out.append((x, y))
    return sorted(out)


def _oracle_negpell(D, N, xhi):
    """x^2 - D y^2 == N, x>0, y>0, x<=xhi."""
    out = []
    y = 1
    while D * y * y + (N if N < 0 else 0) <= xhi * xhi:
        v = N + D * y * y
        if v > 0 and _sq(v):
            x = isqrt(v)
            if 0 < x <= xhi and y > 0:
                out.append((x, y))
        y += 1
        if y > xhi:
            break
    return sorted(out)


CASES = [
    # ---- linear systems (the class the held-out run first broke on) --------
    dict(label="lin-sys-2eq-3var",
         math="Solve[{x + 2 y + 3 z == 10, x - y + z == 2}, {x, y, z}, Integers]",
         vars=["x", "y", "z"], box={"x": (-25, 25), "y": (-25, 25), "z": (-25, 25)},
         sat=lambda d: d["x"] + 2 * d["y"] + 3 * d["z"] == 10 and d["x"] - d["y"] + d["z"] == 2,
         kind="param", source="underdetermined linear system (HNF family)"),
    dict(label="lin-sys-3eq-4var",
         math="Solve[{a + b + c + d == 10, a - b + c - d == 2, a + b - c - d == 0}, {a, b, c, d}, Integers]",
         vars=["a", "b", "c", "d"], box={k: (-14, 14) for k in "abcd"},
         sat=lambda d: (d["a"] + d["b"] + d["c"] + d["d"] == 10
                        and d["a"] - d["b"] + d["c"] - d["d"] == 2
                        and d["a"] + d["b"] - d["c"] - d["d"] == 0),
         kind="param", source="3 eqns / 4 unknowns"),
    dict(label="lin-sys-determined",
         math="Solve[{x + y == 5, x - y == 1}, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-20, 20), "y": (-20, 20)},
         sat=lambda d: d["x"] + d["y"] == 5 and d["x"] - d["y"] == 1,
         kind="finite", source="determined system, unique integer solution"),
    dict(label="lin-sys-noninteger",
         math="Solve[{x + y == 2, x - y == 1}, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-20, 20), "y": (-20, 20)},
         sat=lambda d: d["x"] + d["y"] == 2 and d["x"] - d["y"] == 1,
         kind="empty", source="determined but half-integer -> {} proof"),
    dict(label="lin-sys-int-inconsistent",
         math="Solve[{2 x + 2 y == 3, x - y == 0}, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-40, 40), "y": (-40, 40)},
         sat=lambda d: 2 * d["x"] + 2 * d["y"] == 3 and d["x"] - d["y"] == 0,
         kind="empty", source="parity obstruction -> {} proof"),

    # ---- factorable binary quadratic (Runge) -------------------------------
    dict(label="factorable-cross-4",
         math="Solve[x^2 + x y - 2 y^2 == 4, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-30, 30), "y": (-30, 30)},
         sat=lambda d: d["x"] ** 2 + d["x"] * d["y"] - 2 * d["y"] ** 2 == 4,
         kind="finite", source="(x-y)(x+2y)==4, six points"),
    dict(label="factorable-prod-15-empty",
         math="Solve[(x - y) (x + 2 y) == 15, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-60, 60), "y": (-60, 60)},
         sat=lambda d: (d["x"] - d["y"]) * (d["x"] + 2 * d["y"]) == 15,
         kind="empty", source="mod-3 obstruction -> {} proof"),
    dict(label="factorable-nonunit-7",
         math="Solve[2 x^2 + 3 x y - 2 y^2 == 7, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-40, 40), "y": (-40, 40)},
         sat=lambda d: 2 * d["x"] ** 2 + 3 * d["x"] * d["y"] - 2 * d["y"] ** 2 == 7,
         kind="finite", source="(2x-y)(x+2y)==7, non-unit leading coeffs"),
    dict(label="diff-of-squares-21",
         math="Solve[x^2 - y^2 == 21, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-30, 30), "y": (-30, 30)},
         sat=lambda d: d["x"] ** 2 - d["y"] ** 2 == 21,
         kind="finite", source="difference of squares"),

    # ---- definite binary quadratic / rotated ellipse (Tier-2 D, Delta<0) ---
    dict(label="ellipse-rotated-7",
         math="Solve[x^2 + x y + y^2 == 7, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-6, 6), "y": (-6, 6)},
         sat=lambda d: d["x"] ** 2 + d["x"] * d["y"] + d["y"] ** 2 == 7,
         kind="finite", source="rotated ellipse, cross term, 12 points"),
    dict(label="ellipse-nonunit-24",
         math="Solve[3 x^2 + 2 x y + 3 y^2 == 24, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-6, 6), "y": (-6, 6)},
         sat=lambda d: 3 * d["x"] ** 2 + 2 * d["x"] * d["y"] + 3 * d["y"] ** 2 == 24,
         kind="finite", source="non-unit rotated ellipse"),
    dict(label="ellipse-linear-empty",
         math="Solve[x^2 + x y + y^2 - 3 x == 7, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-12, 12), "y": (-12, 12)},
         sat=lambda d: d["x"] ** 2 + d["x"] * d["y"] + d["y"] ** 2 - 3 * d["x"] == 7,
         kind="empty", source="definite form + linear terms, no solution -> {} proof"),
    dict(label="ellipse-norep-empty",
         math="Solve[5 x^2 + 6 x y + 5 y^2 == 8, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-20, 20), "y": (-20, 20)},
         sat=lambda d: 5 * d["x"] ** 2 + 6 * d["x"] * d["y"] + 5 * d["y"] ** 2 == 8,
         kind="empty", source="8 not represented by the form -> {} proof"),

    # ---- generalised Pell (unbounded, positive orthant) --------------------
    dict(label="genpell-2-7",
         math="Solve[x^2 - 2 y^2 == 7 && x > 0 && y > 0, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (1, 400), "y": (1, 400)},
         sat=lambda d: d["x"] ** 2 - 2 * d["y"] ** 2 == 7 and d["x"] > 0 and d["y"] > 0,
         kind="param", source="two classes"),
    dict(label="genpell-5-4",
         math="Solve[x^2 - 5 y^2 == 4 && x > 0 && y > 0, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (1, 400), "y": (1, 400)},
         sat=lambda d: d["x"] ** 2 - 5 * d["y"] ** 2 == 4 and d["x"] > 0 and d["y"] > 0,
         kind="param", source="three classes"),
    dict(label="genpell-2-5-empty",
         math="Solve[x^2 - 2 y^2 == 5 && x > 0 && y > 0, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (1, 400), "y": (1, 400)},
         sat=lambda d: d["x"] ** 2 - 2 * d["y"] ** 2 == 5 and d["x"] > 0 and d["y"] > 0,
         kind="empty", source="mod-8 obstruction -> {} proof"),
    dict(label="negpell-13-solvable",
         math="Solve[x^2 - 13 y^2 == -1 && x > 0 && y > 0, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (1, 4000), "y": (1, 4000)},
         sat=lambda d: d["x"] ** 2 - 13 * d["y"] ** 2 == -1 and d["x"] > 0 and d["y"] > 0,
         kind="param", oracle=lambda: _oracle_negpell(13, -1, 4000), source="negative Pell, one class"),
    dict(label="negpell-3-empty",
         math="Solve[x^2 - 3 y^2 == -1 && x > 0 && y > 0, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (1, 4000), "y": (1, 4000)},
         sat=lambda d: d["x"] ** 2 - 3 * d["y"] ** 2 == -1 and d["x"] > 0 and d["y"] > 0,
         kind="empty", oracle=lambda: _oracle_negpell(3, -1, 4000), source="even CF period -> {} proof"),

    # ---- bounded nonlinear (the strong existing engine) --------------------
    dict(label="mordell-neg2-boxed",
         math="Solve[y^2 == x^3 - 2 && 0 < x < 1000, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (1, 999), "y": (-4000, 4000)},
         sat=lambda d: d["y"] ** 2 == d["x"] ** 3 - 2 and 0 < d["x"] < 1000,
         kind="finite", oracle=lambda: _oracle_mordell(-2, 1, 999), source="Mordell k=-2 (3,+/-5)"),
    dict(label="thue-cubic-boxed",
         math="Solve[x^3 - 2 y^3 == 1 && Abs[x] < 500 && Abs[y] < 500, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-499, 499), "y": (-499, 499)},
         sat=lambda d: d["x"] ** 3 - 2 * d["y"] ** 3 == 1 and abs(d["x"]) < 500 and abs(d["y"]) < 500,
         kind="finite", oracle=lambda: _oracle_thue_x3_m2y3(500), source="Thue x^3-2y^3=1, bounded"),
    dict(label="thue-cubic-m3-unbounded",
         math="Solve[x^3 - 2 y^3 == 3, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-200, 200), "y": (-200, 200)},
         sat=lambda d: d["x"] ** 3 - 2 * d["y"] ** 3 == 3,
         kind="finite", source="Thue x^3-2y^3=3 unbounded (M2 general-m); complete set {(-5,-4),(1,-1)}"),
    dict(label="thue-cubic-m5-empty",
         math="Solve[x^3 - 2 y^3 == 5, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-200, 200), "y": (-200, 200)},
         sat=lambda d: d["x"] ** 3 - 2 * d["y"] ** 3 == 5,
         kind="empty", source="Thue x^3-2y^3=5 has no norm-5 element -> proven {} (not merely a decline)"),
    dict(label="thue-cyclo-phi5-unbounded",
         math="Solve[x^4 + x^3 y + x^2 y^2 + x y^3 + y^4 == 1, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-20, 20), "y": (-20, 20)},
         sat=lambda d: d["x"] ** 4 + d["x"] ** 3 * d["y"] + d["x"] ** 2 * d["y"] ** 2
                       + d["x"] * d["y"] ** 3 + d["y"] ** 4 == 1,
         kind="finite", source="cyclotomic Phi5 Thue over Q(zeta5) (totally complex, r1=0); 6 solutions"),
    dict(label="thue-totallycomplex-x4y4-17",
         math="Solve[x^4 + y^4 == 17, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-20, 20), "y": (-20, 20)},
         sat=lambda d: d["x"] ** 4 + d["y"] ** 4 == 17,
         kind="finite", source="x^4+y^4=17 over Q(zeta8) (totally complex, general m); 8 solutions"),
    dict(label="thue-totallycomplex-x4y4-3-empty",
         math="Solve[x^4 + y^4 == 3, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-20, 20), "y": (-20, 20)},
         sat=lambda d: d["x"] ** 4 + d["y"] ** 4 == 3,
         kind="empty", source="x^4+y^4=3 (not a sum of two 4th powers) -> proven {} via the |Im| bound"),
    dict(label="thue-rank2-quartic-d10",
         math="Solve[x^4 - 10 y^4 == 1, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-50, 50), "y": (-50, 50)},
         sat=lambda d: d["x"] ** 4 - 10 * d["y"] ** 4 == 1,
         kind="finite", source="Q(10^1/4) rank-2 quartic (M4 Voronoi minima walk); complete set {(-1,0),(1,0)}"),
    dict(label="thue-rank2-quartic-d10-neg-empty",
         math="Solve[x^4 - 10 y^4 == -1, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-50, 50), "y": (-50, 50)},
         sat=lambda d: d["x"] ** 4 - 10 * d["y"] ** 4 == -1,
         kind="empty", source="Q(10^1/4) rank-2 quartic; x^4-10y^4=-1 proven {}"),
    dict(label="two-squares-25",
         math="Solve[x^2 + y^2 == 25, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-6, 6), "y": (-6, 6)},
         sat=lambda d: d["x"] ** 2 + d["y"] ** 2 == 25,
         kind="finite", source="sum of two squares, 12 signed pairs"),
    dict(label="ternary-legendre-empty",
         math=("Solve[3 x^2 + 5 y^2 == 7 z^2 && Abs[x] + Abs[y] + Abs[z] > 0 "
               "&& Abs[x] < 40 && Abs[y] < 40 && Abs[z] < 40, {x, y, z}, Integers]"),
         vars=["x", "y", "z"], box={"x": (-39, 39), "y": (-39, 39), "z": (-39, 39)},
         sat=lambda d: (3 * d["x"] ** 2 + 5 * d["y"] ** 2 == 7 * d["z"] ** 2
                        and abs(d["x"]) + abs(d["y"]) + abs(d["z"]) > 0),
         kind="empty", source="Legendre ternary, no nontrivial solution in box"),

    # ---- genuinely research-grade: a decline is acceptable, a WRONG is not --
    dict(label="bqf-nonsquare-disc-unbounded",
         math="Solve[x^2 + 3 x y + y^2 == 11, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-60, 60), "y": (-60, 60)},
         sat=lambda d: d["x"] ** 2 + 3 * d["x"] * d["y"] + d["y"] ** 2 == 11,
         kind="decline_ok", source="Pell-type conic (delta=5); unbounded -> decline OK"),
    dict(label="mordell-pos-unbounded",
         math="Solve[y^2 == x^3 + 17, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-3, 60), "y": (-500, 500)},
         sat=lambda d: d["y"] ** 2 == d["x"] ** 3 + 17,
         kind="decline_ok", oracle=lambda: _oracle_mordell(17, -3, 60), source="elliptic integral points; decline OK, must not be wrong {}"),
    dict(label="thue-quartic-generalm-unbounded",
         math="Solve[x^4 - 3 y^4 == 13, {x, y}, Integers]",
         vars=["x", "y"], box={"x": (-50, 50), "y": (-50, 50)},
         sat=lambda d: d["x"] ** 4 - 3 * d["y"] ** 4 == 13,
         kind="decline_ok", source="quartic Thue, |m|!=1 (out of the cubic mu-scope) -> decline OK; any answer must be correct"),
]
