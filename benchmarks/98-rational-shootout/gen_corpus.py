#!/usr/bin/env python3.11
"""Generate the Together / Cancel / Apart shootout corpus (deterministic, seeded).

Writes corpus.json: {"seed", "families": [...], "cases": [...]} where every case is

    {"id", "family", "op": "together" | "cancel" | "apart", "size",
     "vars": [...], "var": "x" (apart), "input": "<text>",
     "terms": [...]   (together: the summands, so PARI can pre-parse them),
     "num", "den"     (cancel: the two sides, ditto),
     "options": {} | {"modulus": p} | {"extension": "auto"}}

Text is a neutral syntax every engine reads after light rewriting: integers,
variable names, `+ - * / ^`, parentheses, `I` for the imaginary unit,
`sqrt(2)`, and two kernel tokens: `LG` = log(x) and `EX` = exp(x).  Polynomial
pieces are written expanded.

    python3.11 gen_corpus.py          # writes corpus.json next to this file
"""
import os
import sys

try:
    import sympy  # noqa: F401
except ImportError:                      # re-exec under the Python that has sympy
    os.execv("/usr/local/bin/python3.11", ["/usr/local/bin/python3.11"] + sys.argv)

import itertools
import json
import random

from sympy import I, Poly, ZZ, expand, prod, sqrt, symbols

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20261028
rng = random.Random(SEED)
P61 = 2305843009213693951            # 2^61 - 1 (prime)

X, Y, Z, W, U, V = symbols("x y z w u v")
A, B, C, D, F, G = symbols("a b c d f g")
LG, EX = symbols("LG EX")             # kernel placeholders: log(x), exp(x)

FAMILIES = []
cases = []


def text(e):
    s = str(expand(e)).replace("**", "^")
    return s.replace("sqrt(2)", "sqrt(2)")


def fam(name):
    FAMILIES.append(name)
    return name


def add(family, op, size, gens, options=None, input_text=None, terms=None, num=None,
        den=None, var=None):
    cid = f"{family}-{len([c for c in cases if c['family'] == family]) + 1}"
    rec = {"id": cid, "family": family, "op": op, "size": size,
           "vars": [str(g) for g in gens], "options": options or {}}
    if op == "together":
        rec["terms"] = [f"({text(n)})/({text(d)})" if d != 1 else f"({text(n)})" for n, d in terms]
        rec["input"] = " + ".join(rec["terms"])
    elif op == "cancel":
        rec["num"], rec["den"] = text(num), text(den)
        rec["input"] = f"({rec['num']})/({rec['den']})"
    else:
        rec["var"] = str(var)
        rec["input"] = input_text if input_text is not None else f"({text(num)})/({text(den)})"
    cases.append(rec)


def rcoef(bits):
    return rng.choice([-1, 1]) * rng.randrange(1, 2 ** bits + 1)


def rand_uni(d, bits, x=X):
    cs = [rcoef(bits)] + [rng.randrange(-2 ** bits, 2 ** bits + 1) for _ in range(d - 1)] + [rcoef(bits)]
    return Poly(cs, x, domain=ZZ).as_expr()


def rand_mv(gens, deg, bits, nterms=None):
    mons = [m for m in itertools.product(range(deg + 1), repeat=len(gens)) if sum(m) <= deg]
    top = [m for m in mons if sum(m) == deg]
    pick = mons if nterms is None else rng.sample(mons, min(nterms, len(mons)))
    if not any(sum(m) == deg for m in pick):
        pick = pick + [rng.choice(top)]
    if not any(sum(m) == 0 for m in pick):
        pick = pick + [(0,) * len(gens)]
    return sum(rcoef(bits) * prod(g ** e for g, e in zip(gens, m)) for m in pick)


# ---------------------------------------------------------------- Together
f = fam("tog-lin")
for n in (20, 50, 100, 200):
    add(f, "together", f"n={n}", [X], terms=[(k, X + k) for k in range(1, n + 1)])

f = fam("tog-uni")
for n in (8, 25, 60, 120):
    ts = []
    for _ in range(n):
        q = X ** 2 + rng.randrange(-9, 10) * X + rng.randrange(1, 30)
        ts.append((rng.randrange(-20, 21) or 1, q ** rng.randrange(1, 3)))
    add(f, "together", f"n={n}", [X], terms=ts)

f = fam("tog-multi")
for n in (4, 8, 12, 16):
    ts = [(rand_mv([X, Y, Z], 2, 4), rand_mv([X, Y, Z], 2, 4)) for _ in range(n)]
    add(f, "together", f"n={n}, 3 vars", [X, Y, Z], terms=ts)

f = fam("tog-common")
for n, gd in ((5, 10), (10, 20), (20, 30)):
    Gc = rand_uni(gd, 6)
    ts = [(rand_uni(3, 6), Gc * (X ** 2 + rng.randrange(-5, 6) * X + rng.randrange(1, 20))) for _ in range(n)]
    add(f, "together", f"n={n}, shared deg {gd}", [X], terms=ts)
for n in (4, 8):
    Gc = rand_mv([X, Y], 3, 4)
    ts = [(rand_mv([X, Y], 2, 4), Gc * rand_mv([X, Y], 1, 4)) for _ in range(n)]
    add(f, "together", f"n={n}, 2 vars", [X, Y], terms=ts)

f = fam("tog-telescope")
for n in (25, 100, 300):
    ts = []
    for k in range(1, n + 1):
        ts += [(1, X + k), (-1, X + k + 1)]
    add(f, "together", f"2n={2 * n}", [X], terms=ts)
for n in (10, 40):
    ts = []
    for k in range(1, n + 1):
        ts += [(1, X + k * Y), (-1, X + (k + 1) * Y)]
    add(f, "together", f"2n={2 * n}, 2 vars", [X, Y], terms=ts)

f = fam("tog-nested")
for depth, gens in ((8, [X]), (16, [X]), (32, [X]), (8, [X, Y]), (14, [X, Y])):
    s = "x"
    for k in range(depth):
        s = f"(x + {k + 1}/({s}))" if len(gens) == 1 else f"(x + y/({s}) + {k + 1})"
    cases.append({"id": f"{f}-{len([c for c in cases if c['family'] == f]) + 1}", "family": f,
                  "op": "together", "size": f"depth {depth}, {len(gens)} var",
                  "vars": [str(g) for g in gens], "options": {}, "terms": [f"1/{s}"],
                  "input": f"1/{s}", "nested": True})

f = fam("tog-kernel")
for n in (4, 8, 16):
    ts = []
    for _ in range(n):
        nu = rng.randrange(-5, 6) * LG + rng.randrange(-5, 6) * X + rng.randrange(1, 6)
        de = rng.choice([EX + rng.randrange(1, 5), EX ** 2 - rng.randrange(1, 5) ** 2,
                         EX * X + rng.randrange(1, 5), LG + X + rng.randrange(1, 5)])
        ts.append((nu, de))
    add(f, "together", f"n={n}", [X, LG, EX], terms=ts)

f = fam("tog-gauss")
for n in (4, 8, 16):
    ts = [(rng.randrange(1, 9), X + rng.randrange(-5, 6) + rng.randrange(1, 6) * I) for _ in range(n)]
    ts += [(1, X ** 2 + 1)]
    add(f, "together", f"n={n + 1}", [X], terms=ts)

f = fam("tog-mod")
for p in (7, 1000003, P61):
    add(f, "together", f"n=60, p={'2^61-1' if p == P61 else p}", [X], options={"modulus": p},
        terms=[(k, X + k) for k in range(1, 61)])

# ---------------------------------------------------------------- Cancel
f = fam("can-uni")
for gd in (50, 200, 500, 1000):
    Gc = rand_uni(gd, 8)
    add(f, "cancel", f"gcd deg {gd}", [X], num=rand_uni(20, 8) * Gc, den=rand_uni(20, 8) * Gc)

f = fam("can-uni-big")
for gd in (20, 80, 200):
    Gc = rand_uni(gd, 200)
    add(f, "cancel", f"gcd deg {gd}, 200-bit", [X], num=rand_uni(10, 200) * Gc, den=rand_uni(10, 200) * Gc)

f = fam("can-multi-dense")
for gd in (4, 6, 8):
    Gc = rand_mv([X, Y, Z], gd, 6)
    add(f, "cancel", f"3 vars, gcd deg {gd}", [X, Y, Z],
        num=rand_mv([X, Y, Z], 3, 6) * Gc, den=rand_mv([X, Y, Z], 3, 6) * Gc)

f = fam("can-multi-sparse")
g6 = [X, Y, Z, W, U, V]
for nt in (10, 20, 40):
    Gc = rand_mv(g6, 8, 8, nterms=nt)
    add(f, "cancel", f"6 vars, gcd {nt} terms", g6,
        num=rand_mv(g6, 6, 8, nterms=nt // 2) * Gc, den=rand_mv(g6, 6, 8, nterms=nt // 2) * Gc)

f = fam("can-coprime")
for d in (100, 400, 1000):
    add(f, "cancel", f"deg {d}", [X], num=rand_uni(d, 8), den=rand_uni(d, 8))
for d in (5, 8):
    add(f, "cancel", f"3 vars deg {d}", [X, Y, Z], num=rand_mv([X, Y, Z], d, 6), den=rand_mv([X, Y, Z], d, 6))

f = fam("can-small")
for nu, de, gens in ((X ** 2 - 1, X - 1, [X]), (X ** 3 - Y ** 3, X - Y, [X, Y]),
                     (X ** 4 - 1, X ** 6 - 1, [X]), (A ** 2 - B ** 2, A + B, [A, B]),
                     ((X + 1) ** 3 * (X - 2), (X + 1) * (X - 2) ** 2 * (X + 3), [X]),
                     (X ** 2 * Y - Y, X * Y + Y, [X, Y])):
    add(f, "cancel", "textbook", gens, num=nu, den=de)

f = fam("can-gauss")
for d in (5, 20, 60):
    Gc = prod(X - (rng.randrange(-9, 10) + rng.randrange(1, 9) * I) for _ in range(d))
    add(f, "cancel", f"gcd deg {d}", [X], num=rand_uni(5, 6) * Gc, den=rand_uni(5, 6) * Gc)

f = fam("can-ext")
for d in (5, 20, 60):
    add(f, "cancel", f"cofactor deg {d}", [X], options={"extension": "auto"},
        num=(X ** 2 - 2) * rand_uni(d, 6), den=(X - sqrt(2)) * rand_uni(d, 6))

f = fam("can-mod")
for p, gd in ((7, 200), (1000003, 200), (P61, 200)):
    Gc = rand_uni(gd, 8)
    add(f, "cancel", f"gcd deg {gd}, p={'2^61-1' if p == P61 else p}", [X], options={"modulus": p},
        num=rand_uni(20, 8) * Gc, den=rand_uni(20, 8) * Gc)

# ---------------------------------------------------------------- Apart
f = fam("ap-lin")
for n in (10, 30, 60, 100):
    add(f, "apart", f"{n} linear", [X], var=X, num=rand_uni(n - 1, 6), den=prod(X + k for k in range(1, n + 1)))

f = fam("ap-mult")
for a, b, c in ((3, 3, 3), (10, 8, 5), (30, 20, 10), (60, 40, 20)):
    add(f, "apart", f"mult {a},{b},{c}", [X], var=X, num=rand_uni(5, 6),
        den=(X + 1) ** a * (X - 2) ** b * (2 * X + 3) ** c)

f = fam("ap-irr")
qs = [X ** 2 + 1, X ** 2 + X + 1, X ** 3 + 2, 3 * X ** 2 + 5, X ** 4 + X + 1, X ** 2 - 3]
for k, m in ((2, 1), (4, 2), (6, 3)):
    add(f, "apart", f"{k} irreducibles, mult {m}", [X], var=X, num=rand_uni(4, 6),
        den=prod(q ** m for q in qs[:k]))

f = fam("ap-improper")
for n in (5, 20, 50):
    add(f, "apart", f"{n} linear, num deg {n + 10}", [X], var=X, num=rand_uni(n + 10, 6),
        den=prod(X - k for k in range(1, n + 1)))

f = fam("ap-param")
ps = [A, B, C, D, F, G]
for n in (2, 3, 4, 5, 6):
    add(f, "apart", f"{n} params", [X] + ps[:n], var=X, num=1, den=prod(X + p for p in ps[:n]))
add(f, "apart", "(x+a)^2 (x^2+b)", [X, A, B], var=X, num=X + 1, den=(X + A) ** 2 * (X ** 2 + B))

f = fam("ap-param-content")
for n in (4, 10, 20):
    add(f, "apart", f"{n} linear, param in num/content", [X, A], var=X,
        num=rand_uni(n - 1, 4) + A * X, den=A * prod(X - k for k in range(1, n + 1)))

f = fam("ap-gauss")
for n in (3, 6, 10):
    add(f, "apart", f"{n} gaussian linear", [X], var=X, num=1,
        den=prod(X - (k + (k % 3 + 1) * I) for k in range(1, n + 1)))

f = fam("ap-mod")
for p in (7, 1000003, P61):
    add(f, "apart", f"6 linear, p={'2^61-1' if p == P61 else p}", [X], var=X, options={"modulus": p},
        num=rand_uni(4, 6), den=prod(X + k for k in range(1, 7)))

json.dump({"seed": SEED, "families": FAMILIES, "cases": cases},
          open(os.path.join(HERE, "corpus.json"), "w"), indent=1)
print(f"{len(cases)} cases in {len(FAMILIES)} families")
