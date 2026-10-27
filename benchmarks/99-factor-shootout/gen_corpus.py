#!/usr/bin/env python3.11
"""Generate the polynomial Factor shootout corpus (deterministic, seeded).

Writes corpus.json: {"seed": ..., "cases": [...]} where every case is

    {"id", "family", "size", "vars": [...], "input": "<expanded text>",
     "options": {} | {"modulus": p} | {"gaussian": true} | {"extension": [2, 3]}}

The input is always given EXPANDED (the product is multiplied out before it is
written), in a neutral syntax every engine reads unchanged: integers, variable
names, `+ - * ^`, parentheses, and for a rational function `(num)/(den)`.
`extension: [2, 3]` means the field Q(Sqrt[2], Sqrt[3]).

Families are listed in FAMILIES (order = report order).  Sizes are chosen so
each case finishes well under ~30 s in the fastest engine.

The Swinnerton-Dyer polynomials are computed with PARI (`gp`) by iterated
resultants, because sympy needs minutes for degree 64.

    python3.11 gen_corpus.py          # writes corpus.json next to this file
"""
import os
import sys

try:
    import sympy  # noqa: F401
except ImportError:                      # re-exec under the Python that has sympy
    os.execv("/usr/local/bin/python3.11", ["/usr/local/bin/python3.11"] + sys.argv)

import json
import random
import subprocess

from sympy import Matrix, Poly, ZZ, cancel, fraction, symbols, together

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20261027
rng = random.Random(SEED)
P60 = 1152921504606846883            # a 60-bit prime (2^60 - 93)

X, Y, Z, W, U, V = symbols("x y z w u v")
VARS6 = [X, Y, Z, W, U, V]

cases = []


def text(e):
    """Neutral text of a polynomial (or rational) expression."""
    return str(e).replace("**", "^")


def add(family, size, poly, gens, options=None, den=None):
    """Record one case.  `poly` / `den` are sympy Poly or expressions."""
    num = poly.as_expr() if isinstance(poly, Poly) else poly
    if den is not None:
        d = den.as_expr() if isinstance(den, Poly) else den
        s = f"({text(num.expand())})/({text(d.expand())})"
    else:
        s = text(num.expand())
    cid = f"{family}-{len([c for c in cases if c['family'] == family]) + 1}"
    cases.append({"id": cid, "family": family, "size": size,
                  "vars": [str(g) for g in gens], "input": s,
                  "options": options or {}})


def rcoef(bits):
    return rng.choice([-1, 1]) * rng.randrange(1, 2 ** bits + 1)


def rand_uni(d, bits, x=X):
    """Dense univariate of exact degree d with nonzero constant term."""
    cs = [rcoef(bits)] + [rng.randrange(-2 ** bits, 2 ** bits + 1) for _ in range(d - 1)] + [rcoef(bits)]
    return Poly(cs, x, domain=ZZ)


def rand_mv(gens, deg, bits, density=1.0, nterms=None):
    """Random multivariate of total degree <= deg (exact total degree deg)."""
    import itertools
    mons = [m for m in itertools.product(range(deg + 1), repeat=len(gens)) if sum(m) <= deg]
    if nterms is not None:
        top = [m for m in mons if sum(m) == deg]
        pick = rng.sample(mons, min(nterms - 1, len(mons))) + [rng.choice(top)]
    else:
        pick = [m for m in mons if rng.random() < density]
        pick.append(tuple([deg] + [0] * (len(gens) - 1)))
    pick.append(tuple([0] * len(gens)))                 # constant term
    terms = {m: rcoef(bits) for m in pick}
    return Poly.from_dict(terms, *gens, domain=ZZ)


def prod(ps):
    out = ps[0]
    for p in ps[1:]:
        out = out * p
    return out


def swinnerton_dyer(n):
    """Swinnerton-Dyer polynomial for the first n primes, via PARI."""
    primes = [2, 3, 5, 7, 11, 13, 17, 19][:n]
    script = ("f=x; foreach(%s, p, f=polresultant(subst(f,x,x-t), t^2-p, t)); print(f);\nquit;\n"
              % str(primes))
    out = subprocess.run(["gp", "-q", "-D", "parisizemax=1000000000"], input=script,
                         capture_output=True, text=True, check=True).stdout
    from sympy import sympify
    return Poly(sympify(out.strip().replace("^", "**"), locals={"x": X}), X)


# --------------------------------------------------------------------------
# Univariate over Z
# --------------------------------------------------------------------------
for n in (120, 360, 840, 1260, 2310):
    add("cyclotomic", f"x^{n}-1", Poly(X ** n - 1, X), [X])
add("cyclotomic", "x^1155+1", Poly(X ** 1155 + 1, X), [X])

for n in (4, 5, 6):
    sd = swinnerton_dyer(n)
    add("swinnerton-dyer", f"deg {sd.degree()}", sd, [X])
sd4 = swinnerton_dyer(4)
add("swinnerton-dyer", "SD16(x)*SD16(x+1)", sd4 * Poly(sd4.as_expr().subs(X, X + 1), X), [X])

for k, d, b in ((4, 8, 4), (6, 15, 8), (8, 25, 8), (5, 60, 10), (10, 40, 16)):
    add("uni-dense", f"{k}x deg {d}, {b}-bit", prod([rand_uni(d, b) for _ in range(k)]), [X])

for k, d, b in ((2, 10, 100), (3, 15, 300), (3, 30, 600), (4, 40, 1000)):
    add("uni-bigcoef", f"{k}x deg {d}, {b}-bit", prod([rand_uni(d, b) for _ in range(k)]), [X])

for spec in (((3, 10), (5, 6)), ((2, 20), (4, 10), (7, 5)), ((5, 15), (9, 8), (12, 4), (3, 25))):
    ps = [rand_uni(d, 6) ** m for d, m in spec]
    add("uni-multiplicity", " ".join(f"d{d}^{m}" for d, m in spec), prod(ps), [X])


def sparse_uni(terms, maxdeg):
    exps = sorted(rng.sample(range(1, maxdeg), terms - 2)) + [maxdeg]
    return Poly(sum(rcoef(4) * X ** e for e in exps) + rcoef(4), X)


for k, t, D in ((2, 4, 300), (3, 3, 600), (2, 5, 1500), (3, 4, 2500)):
    add("uni-sparse", f"{k}x {t}-term deg {D}", prod([sparse_uni(t, D) for _ in range(k)]), [X])

# --------------------------------------------------------------------------
# Multivariate over Q
# --------------------------------------------------------------------------
for nv, k, d, b in ((2, 3, 5, 8), (2, 4, 10, 8), (3, 3, 4, 6), (3, 4, 6, 6), (4, 3, 3, 4)):
    g = VARS6[:nv]
    add("mv-dense", f"{nv}v {k}x deg {d}", prod([rand_mv(g, d, b) for _ in range(k)]), g)

for nv, k, t, d in ((4, 3, 5, 6), (5, 3, 6, 8), (6, 3, 6, 10), (6, 4, 8, 12)):
    g = VARS6[:nv]
    add("mv-sparse", f"{nv}v {k}x {t}-term deg {d}",
        prod([rand_mv(g, d, 5, nterms=t) for _ in range(k)]), g)


def wang_factor(gens, d, lcdeg):
    """x^d * lc(y, z, ...) + lower-degree terms: non-monic in the main variable."""
    x, rest = gens[0], gens[1:]
    lc = rand_mv(rest, lcdeg, 4, nterms=3).as_expr()
    low = rand_mv(gens, d - 1, 5, nterms=5).as_expr()
    return Poly(lc * x ** d + low, *gens)


for nv, k, d, l in ((3, 2, 3, 2), (3, 3, 4, 2), (4, 3, 4, 3), (5, 3, 3, 3)):
    g = VARS6[:nv]
    add("mv-nonmonic", f"{nv}v {k}x deg {d}, lc deg {l}", prod([wang_factor(g, d, l) for _ in range(k)]), g)

for nv, d in ((2, 20), (3, 8), (4, 6), (5, 4)):
    g = VARS6[:nv]
    add("mv-irreducible", f"{nv}v deg {d}", rand_mv(g, d, 8, density=0.7), g)

for n, sign in ((6, -1), (10, 1), (15, 1), (20, 1)):
    add("mv-power", f"(x+y+z)^{n}{'+' if sign > 0 else '-'}1", Poly((X + Y + Z) ** n + sign, X, Y, Z), [X, Y, Z])
add("mv-power", "(x+y+z+w)^8-1", Poly((X + Y + Z + W) ** 8 - 1, X, Y, Z, W), [X, Y, Z, W])

for n in (4, 5, 6, 7):
    g = list(symbols(" ".join(f"x{i}" for i in range(1, n + 1))))
    vdm = Matrix(n, n, lambda i, j: g[j] ** i)
    add("mv-determinant", f"Vandermonde {n}", Poly(vdm.det(method="berkowitz").expand(), *g), g)
for n in (3, 4):
    xs = list(symbols(" ".join(f"x{i}" for i in range(1, n + 1))))
    ys = list(symbols(" ".join(f"y{i}" for i in range(1, n + 1))))
    num = prod([Poly((xs[j] - xs[i]) * (ys[j] - ys[i]), *(xs + ys)) for i in range(n) for j in range(i + 1, n)])
    add("mv-determinant", f"Cauchy numerator {n}", num, xs + ys)

# --------------------------------------------------------------------------
# Rational functions
# --------------------------------------------------------------------------
for d, k in ((10, 2), (25, 3), (40, 3)):
    a, b, c = (prod([rand_uni(d, 6) for _ in range(k)]) for _ in range(3))
    add("rational", f"uni deg {d * k * 2}", a * b, [X], den=a * c)
for nv, d in ((2, 4), (3, 3)):
    g = VARS6[:nv]
    a, b, c = (rand_mv(g, d, 5) for _ in range(3))
    add("rational", f"{nv}v deg {2 * d}", a * b, g, den=a * c * c)
xs = list(symbols("x1 x2 x3"))
ys = list(symbols("y1 y2 y3"))
cdet = together(Matrix(3, 3, lambda i, j: 1 / (xs[i] + ys[j])).det(method="berkowitz"))
cn, cd = fraction(cancel(cdet))
add("rational", "Cauchy det 3", cn, xs + ys, den=cd)

# --------------------------------------------------------------------------
# Modulus p (univariate and multivariate)
# --------------------------------------------------------------------------
for p, d, k in ((2, 1000, 3), (3, 500, 3), (10007, 300, 3), (P60, 150, 3)):
    ps = [rand_uni(d // k, 30) for _ in range(k)] + [rand_uni(d // 4, 30)]
    add("modulus-uni", f"p={p if p < 10**6 else '2^60-93'} deg {d + d // 4}",
        prod(ps), [X], {"modulus": p})
for p, nv, d, k in ((2, 2, 6, 3), (3, 3, 4, 3), (10007, 2, 8, 3), (P60, 3, 3, 3)):
    g = VARS6[:nv]
    add("modulus-mv", f"p={p if p < 10**6 else '2^60-93'} {nv}v {k}x deg {d}",
        prod([rand_mv(g, d, 20) for _ in range(k)]), g, {"modulus": p})

# --------------------------------------------------------------------------
# GaussianIntegers -> True
# --------------------------------------------------------------------------
add("gaussian", "x^24-1", Poly(X ** 24 - 1, X), [X], {"gaussian": True})
add("gaussian", "x^120-1", Poly(X ** 120 - 1, X), [X], {"gaussian": True})
add("gaussian", "(x^2+k) products",
    prod([Poly(X ** 2 + rng.randrange(1, 50) ** 2, X) for _ in range(4)] + [rand_uni(10, 6)]),
    [X], {"gaussian": True})
add("gaussian", "x^2+y^2 products",
    Poly((X ** 2 + Y ** 2) * (X ** 2 + 4 * Y ** 2) * (X + Y + 1), X, Y), [X, Y], {"gaussian": True})
add("gaussian", "3v sums of squares",
    Poly((X ** 2 + Y ** 2) * (Z ** 2 + 9 * X ** 2) * (Y ** 2 + Z ** 2 + X), X, Y, Z), [X, Y, Z],
    {"gaussian": True})

# --------------------------------------------------------------------------
# Extension -> Sqrt[2] / {Sqrt[2], Sqrt[3]}
# --------------------------------------------------------------------------
add("extension", "Sqrt[2] uni deg 16",
    prod([Poly(X ** 2 - 2, X), Poly(X ** 4 - 8, X), rand_uni(10, 6)]), [X], {"extension": [2]})
add("extension", "Sqrt[2] uni deg 46",
    prod([Poly(X ** 4 - 2 * X ** 2 - 1, X), Poly(X ** 2 - 18, X), rand_uni(20, 8), rand_uni(20, 8)]),
    [X], {"extension": [2]})
add("extension", "{Sqrt[2],Sqrt[3]} uni deg 16",
    prod([Poly(X ** 4 - 10 * X ** 2 + 1, X), Poly((X + 1) ** 4 - 10 * (X + 1) ** 2 + 1, X),
          rand_uni(8, 6)]), [X], {"extension": [2, 3]})
add("extension", "{Sqrt[2],Sqrt[3]} uni deg 36",
    prod([Poly(X ** 4 - 10 * X ** 2 + 1, X), Poly(X ** 2 - 6, X), Poly(X ** 2 - 12, X),
          rand_uni(14, 8), rand_uni(14, 8)]), [X], {"extension": [2, 3]})
add("extension", "Sqrt[2] 2v",
    Poly((X ** 2 - 2 * Y ** 2) * (X ** 2 - 8 * Y ** 2 + X) * (X * Y + 3), X, Y), [X, Y], {"extension": [2]})
add("extension", "{Sqrt[2],Sqrt[3]} 3v",
    Poly((X ** 2 - 2 * Y ** 2) * (X ** 2 - 3 * Z ** 2) * (Y ** 2 - 6 * Z ** 2 + X) * (X + Y + Z), X, Y, Z),
    [X, Y, Z], {"extension": [2, 3]})

FAMILIES = []
for c in cases:
    if c["family"] not in FAMILIES:
        FAMILIES.append(c["family"])

with open(os.path.join(HERE, "corpus.json"), "w") as f:
    json.dump({"seed": SEED, "families": FAMILIES, "cases": cases}, f, indent=1)
print(f"{len(cases)} cases in {len(FAMILIES)} families -> corpus.json")
for fam in FAMILIES:
    cs = [c for c in cases if c["family"] == fam]
    print(f"  {fam:18s} {len(cs):2d}  max input {max(len(c['input']) for c in cs):8d} chars")
