#!/usr/bin/env python3.11
"""gen_corpus.py -- builds corpus.tsv for benchmark 266 (rational integration).

Deterministic (fixed seeds): re-running it reproduces corpus.tsv byte for byte.
Every case is a rational function of x.  Columns:

    id  family  integrand(Wolfram InputForm)  sympy(str, or "-")  note

The sympy column is what tools/... feeds to sympy's ratint; "-" marks a case sympy
cannot take as written (a Root object coefficient).

Families (see README.md for what each one stresses):
    classical  bronstein  sqfree  nonsqfree  multiplicity  rootsum  rioboo
    param_generic  param_special  alg_qi  alg_sqrt  alg_nested  alg_root
    largecoef  stress  czichowski  logtoreal
"""
import os
import random
import sys

from sympy import (I, Poly, Rational, Symbol, diff, expand, factor_list, sqrt,
                   symbols, together, cancel, prod)
from sympy.printing.mathematica import mathematica_code as mc

HERE = os.path.dirname(os.path.abspath(__file__))
x = Symbol("x")
a, b, c = symbols("a b c")

CASES = []


def add(fam, f, note="", sym=None):
    """Record one case.  f is a sympy expression, or a (wolfram, sympy) pair of
    strings when sympy cannot represent the integrand (Root objects)."""
    if isinstance(f, tuple):
        wl, sy = f
    else:
        wl = mc(f)
        sy = str(f) if sym is None else sym
    idx = sum(1 for c_ in CASES if c_[1] == fam) + 1
    CASES.append((f"{fam}-{idx:03d}", fam, wl, sy, note))


def rpoly(rng, deg, lo=-5, hi=5, monic=False):
    while True:
        cs = [rng.randint(lo, hi) for _ in range(deg + 1)]
        if monic:
            cs[-1] = 1
        if cs[-1] != 0:
            return sum(cf * x**i for i, cf in enumerate(cs))


def irreducible(p):
    fl = factor_list(p, x)[1]
    return len(fl) == 1 and fl[0][1] == 1 and Poly(fl[0][0], x).degree() == Poly(p, x).degree()


def rirred(rng, deg, lo=-4, hi=4):
    while True:
        p = rpoly(rng, deg, lo, hi, monic=rng.random() < 0.6)
        if irreducible(p):
            return p


def rnum(rng, deg, lo=-6, hi=6):
    if deg < 0:
        return 1
    return rpoly(rng, deg, lo, hi)


# ---------------------------------------------------------------- classical
def classical():
    T = [
        1/(x**2 + 1), 1/(x**2 - 1), x/(x**2 + 1), 1/(x**3 + 1), 1/(x**3 - 1),
        1/(x**4 + 1), 1/(x**4 - 1), x**2/(x**4 + 1), 1/(x**6 + 1), 1/(x*(x + 1)),
        1/(x**2*(x + 1)), (x + 1)/(x**2 + 2*x + 5), (3*x + 2)/(x**2 - 4*x + 13),
        1/(x**2 + x + 1), (x**3 + 1)/(x**2 - 1), x**4/(x**2 + 1), 1/((x - 1)*(x - 2)*(x - 3)),
        (2*x + 3)/((x - 1)**2*(x + 2)), 1/(x**2 + 1)**2, x/(x**2 + 1)**2,
        1/(x*(x**2 + 1)), (x**2 + 2)/(x*(x**2 + 1)**2), 1/(x**4 + x**2 + 1),
        (5*x**2 + 3*x - 2)/(x**3 + 2*x**2), (x**3 - 4*x - 10)/(x**2 - x - 6),
        1/(x**3 + x), (x**2 - x + 6)/(x**3 + 3*x), 1/(x**2 + 4*x + 8),
        (4*x**2 - 3*x + 2)/(4*x**2 - 4*x + 3), (10)/((x - 1)*(x**2 + 9)),
        x**3/(x**2 + 1)**3, (2*x**2 - x + 4)/(x**3 + 4*x), 1/(x**4 - 16),
        (x**5 + x - 1)/(x**3 + 1), 1/(1 + x)**5, x**2/(1 - x)**7,
        (x - 3)/(x**2 + 2*x + 4)**2, 1/(x**5 - x), (x**4 + 1)/(x**5 + 4*x**3),
        1/(x**3 - 2), 1/(x**2 - 2), (1 + x**2)/(1 + x**4), (1 - x**2)/(1 + x**4),
        (x**2)/(x**6 - 1), 1/(x*(x**5 + 1)),
    ]
    for f in T:
        add("classical", f)


# ---------------------------------------------------------------- bronstein
def bronstein():
    B = [
        (1/(x**3 + x), "Ex. 2.1.1"),
        ((x**7 - 24*x**4 - 4*x**2 + 8*x - 8)/(x**8 + 6*x**6 + 12*x**4 + 8*x**2),
         "Ex. 2.2.1-2.3.1 (Hermite / Mack / Horowitz)"),
        ((x**4 - 3*x**2 + 6)/(x**6 - 5*x**4 + 5*x**2 + 4), "Ex. 2.4.1, 2.5.1, 2.6.1, 2.8 (Rioboo)"),
        (36/(x**5 - 2*x**4 - 2*x**3 + 4*x**2 + x - 2), "Ex. 2.5.2, 2.7.2"),
        ((36*x**6 + 126*x**5 + 183*x**4 + Rational(13807, 6)*x**3 - 407*x**2 - Rational(3242, 5)*x
          + Rational(3044, 15)) / ((x**2 + Rational(7, 6)*x + Rational(1, 3))**2 * (x - Rational(2, 5))**3),
         "Hermite example, Bronstein-Salvy partial fractions"),
        ((x**3 + x)/(x**4 - 3*x**2 + 1)/1, "real roots, irrational"),
        ((x**4 - 3*x**2 + 6)/(x**6 - 5*x**4 + 5*x**2 + 4)**2, "Ex. 2.4.1 squared"),
        (1/(x**2 + 1)**4, "Hermite, multiplicity 4"),
        ((7*x**13 + 10*x**8 + 4*x**7 - 7*x**6 - 4*x**3 - 4*x**2 + 3*x + 3)
         / (x**14 - 2*x**8 - 2*x**7 - 2*x**4 - 4*x**3 - x**2 + 2*x + 1), "Lazard-Rioboo / Czichowski stress"),
        ((x**3 - 2*x)/(x**4 - 4*x**2 + 2), "biquadratic real"),
        ((2*x)/(x**4 + 2*x**2 + 5), "ArcTan of a quadratic"),
        ((x**2 + 1)/(x**4 - x**2 + 1), "palindromic"),
        (1/(x**4 + 4), "Sophie Germain"),
        ((5*x**4 + 1)/(x**5 + x + 1)**2, "derivative-divides"),
    ]
    for f, n in B:
        add("bronstein", f, n)


# ---------------------------------------------------------------- squarefree high degree
def sqfree(rng):
    for n in range(5, 17):
        add("sqfree", 1/(x**n - 1), "cyclotomic product")
        add("sqfree", 1/(x**n + 1), "cyclotomic product")
    for k in range(16):
        nf = rng.randint(2, 4)
        den = prod([rirred(rng, rng.randint(2, 5)) for _ in range(nf)])
        dd = Poly(den, x).degree()
        add("sqfree", expand(rnum(rng, dd - 1)) / expand(den), f"deg {dd}")


# ---------------------------------------------------------------- non-squarefree
def nonsqfree(rng):
    for k in range(30):
        nf = rng.randint(1, 3)
        den = prod([rirred(rng, rng.randint(1, 3)) ** rng.randint(2, 4) for _ in range(nf)])
        dd = Poly(den, x).degree()
        f = expand(rnum(rng, dd - 1)) / expand(den)
        add("nonsqfree", f, f"deg {dd}")


# ---------------------------------------------------------------- multiplicity
def multiplicity():
    for n in range(2, 6):
        for k in range(2, 7):
            add("multiplicity", 1/(x**n + 1)**k, f"1/(x^{n}+1)^{k}")
    for k in range(7, 13):
        add("multiplicity", 1/(x**2 + 1)**k)
    for k in range(2, 6):
        add("multiplicity", 1/(x**3 + x + 1)**k)
    for m, k in [(1, 3), (2, 4), (3, 5), (5, 6), (7, 4), (4, 8)]:
        add("multiplicity", x**m/(x**2 + x + 1)**k)


# ---------------------------------------------------------------- RootSum
def rootsum(rng):
    for n in range(3, 14):
        p = x**n + x + 1
        if irreducible(p):
            add("rootsum", 1/p, f"x^{n}+x+1 irreducible")
        p = x**n - x - 1
        if irreducible(p):
            add("rootsum", x/p, f"x^{n}-x-1 irreducible")
    for k in range(14):
        d = rirred(rng, rng.randint(5, 9), -3, 3)
        dd = Poly(d, x).degree()
        add("rootsum", rnum(rng, dd - 1, -3, 3) / d, f"irreducible deg {dd}")


# ---------------------------------------------------------------- Rioboo real forms
def rioboo(rng):
    """d/dx ArcTan[P/Q] = (P'Q - PQ')/(P^2+Q^2): an integrand with a continuous
    real antiderivative whose naive log->ArcTan conversion jumps at every real
    zero of Q.  Plus the classical Rioboo examples."""
    add("rioboo", (x**4 - 3*x**2 + 6)/(x**6 - 5*x**4 + 5*x**2 + 4), "Bronstein 2.8")
    seen = set()
    while sum(1 for c_ in CASES if c_[1] == "rioboo") < 24:
        P = rpoly(rng, rng.randint(1, 4), -3, 3)
        Q = rpoly(rng, rng.randint(1, 4), -3, 3)
        f = cancel(together((diff(P, x)*Q - P*diff(Q, x)) / (P**2 + Q**2)))
        if f == 0 or f.free_symbols != {x}:
            continue
        key = str(f)
        if key in seen:
            continue
        seen.add(key)
        add("rioboo", f, f"D[ArcTan[({P})/({Q})]]")


# ---------------------------------------------------------------- parameters
def param_generic():
    P = [
        1/(x**2 + a), 1/(x**2 - a), 1/(x**2 + a**2), 1/(x**2 - a**2), x/(x**2 + a),
        1/(x**2 + a*x + b), (x + c)/(x**2 + a*x + b), 1/((x - a)*(x - b)), 1/((x - a)*(x - b)*(x - c)),
        1/(x**3 + a**3), 1/(x**3 - a), 1/(x**4 + a**4), x/(x**4 + a**4), 1/(x**4 - a**4),
        1/(x**2 + a)**2, 1/(x**2 + a)**3, x**2/(x**2 + a**2)**2, 1/(a*x**2 + b), 1/(a*x**2 + b*x + c),
        1/(x*(x**2 + a)), 1/((x + a)**2*(x + b)), 1/(x**4 + a*x**2 + b), (a*x + b)/(x**2 + 1),
        1/(x**5 + a), 1/(x**3 + a*x + b), 1/((x**2 + a)*(x**2 + b)), 1/(x - a)**5,
        (x - a)/(x**2 - 2*a*x + a**2 + b**2), 1/(a**2*x**2 - 1), 1/(x**2 + 2*a*x + 1),
        (2*a*x)/(x**4 + a**2), 1/(1 + a**2*x**2)**2, 1/(x**6 + a**6),
    ]
    for f in P:
        add("param_generic", f)


def param_special():
    """A generic family at the parameter values where its generic answer
    degenerates (a zero discriminant, a coinciding root, a vanishing leading
    coefficient)."""
    S = [
        1/(x**2 + 2*a*x + a**2), 1/(x**2 + a*x + a**2/4), 1/((x - a)*(x - a)), 1/((x - a)*(x + a)),
        1/(x**2 + a*x + b).subs({a: 2, b: 1}), 1/(x**2 + a*x + b).subs({a: 0, b: 0}),
        1/(x**2 + a).subs(a, 0), 1/(x**3 + a*x + b).subs({a: -3, b: 2}),
        1/(x**3 + a*x + b).subs({a: 0, b: 0}), 1/((x - a)*(x - b)*(x - c)).subs({c: a}),
        1/(x**4 + a*x**2 + b).subs({a: 2, b: 1}), 1/(x**4 + a*x**2 + a**2/4),
        1/((x**2 + a)*(x**2 + a)), 1/((x**2 + a)*(x**2 + 4*a)), 1/(x**2 + (a - 1)**2),
        1/(x**2 - (a**2 - 1)), 1/(a**2*x**2 + 2*a*x + 1), (x - a)/(x**2 - 2*a*x + a**2),
        x/(x**4 + 2*a**2*x**2 + a**4), 1/(x**3 - 3*a*x**2 + 3*a**2*x - a**3),
    ]
    for f in S:
        add("param_special", f)


# ---------------------------------------------------------------- algebraic coefficients
def alg_qi(rng):
    F = [1/(x**2 + I), 1/(x - I), 1/(x**2 + I*x + 1), 1/(x**3 + I), x/(x**2 + 2*I*x - 2),
         1/((x - I)**2*(x + 1)), (x + I)/(x**2 + 1), 1/(x**4 + I), (1 + I)/(x**2 - I),
         1/(x**2 + 1 + I), 1/((x - I)**3), (2*x + I)/(x**2 + I*x + 3), 1/(x**3 + I*x + 1),
         1/((x**2 + I)**2)]
    for f in F:
        add("alg_qi", f)
    for k in range(6):
        d = expand((x - rng.randint(-3, 3) - rng.randint(-3, 3)*I) * (x**2 + rng.randint(1, 4)*I*x + rng.randint(1, 5)))
        add("alg_qi", rnum(rng, 1, -3, 3) / d)


def alg_sqrt(rng):
    F = [1/(x**2 - sqrt(2)), 1/(x**2 + sqrt(2)), 1/(x - sqrt(3)), 1/(x**2 - 2*sqrt(2)*x + 1),
         1/(x**2 + sqrt(2)*x + 1), x/(x**4 + sqrt(2)), 1/(x**3 - sqrt(2)), 1/((x - sqrt(2))**2*(x + sqrt(3))),
         1/(x**2 + sqrt(5)*x + 2), (sqrt(2)*x + 1)/(x**2 + 2), 1/(x**4 - sqrt(3)*x**2 + 1),
         1/(x**2 - sqrt(2))**2, 1/(x**2 + sqrt(7)), 1/((x**2 - 2)*(x - sqrt(2))),
         sqrt(3)/(x**2 + x + 1), 1/(x**2 + x + sqrt(2))]
    for f in F:
        add("alg_sqrt", f)
    for k in range(4):
        n = rng.choice([2, 3, 5, 6, 7])
        d = expand(x**3 + rng.randint(-2, 2)*sqrt(n)*x + rng.randint(1, 3))
        add("alg_sqrt", 1/d, f"Q(Sqrt[{n}]) cubic")


def alg_nested():
    F = [1/(x**2 - sqrt(2 + sqrt(2))), 1/(x**2 + sqrt(2 + sqrt(3))), 1/(x - sqrt(2) - sqrt(3)),
         1/(x**2 - (sqrt(2) + sqrt(3))), 1/((x - sqrt(2))*(x - sqrt(3))), 1/(x**2 + sqrt(2)*x + sqrt(3)),
         1/(x**4 - 2*sqrt(2)*x**2 - 1), 1/(x**2 - 2**Rational(1, 3)), 1/(x**3 - 2**Rational(2, 3)),
         1/(x**2 + (1 + sqrt(5))/2), x/(x**4 - (1 + sqrt(2))), 1/(x**2 + I*sqrt(2))]
    for f in F:
        add("alg_nested", f)


def alg_root():
    R = [
        ("1/(x - Root[#^3 - 2 &, 1])", "-"),
        ("1/(x^2 - Root[#^3 - #1 - 1 &, 1])", "-"),
        ("1/(x^2 + Root[#^3 + #1 + 1 &, 1] x + 1)", "-"),
        ("1/(x^3 - Root[#^5 - #1 - 1 &, 1])", "-"),
        ("Root[#^3 - 2 &, 1]/(x^2 + 1)", "-"),
        ("1/((x - Root[#^3 - 2 &, 1]) (x - 1))", "-"),
        ("1/(x^2 + Root[#^4 - 2 &, 1])", "-"),
        ("x/(x^4 - Root[#^3 - 3 # + 1 &, 1])", "-"),
        ("1/(x^2 - Root[#^3 - 2 &, 1] x + 1)^2", "-"),
        ("1/(x^3 + Root[837/4 + #1^3 &, 1])", "-"),
        ("With[{r = Root[837/4 + #1^3 &, 1]}, (864 - 2 x^4 r + 32 r^3)/(1109376 - 9 x^6 - 156 x^4 r - 432 x^2 r^2 + 4800 r^3)]",
         "-"),
        ("With[{r = Root[837/4 + #1^3 &, 1]}, x/(x^4 - 9 r x^2 + r^2)]", "-"),
    ]
    notes = [""] * 10 + ["hang repro (Method -> BronsteinRational, > 60 s)", ""]
    for (wl, sy), n in zip(R, notes):
        add("alg_root", (wl, sy), n)


# ---------------------------------------------------------------- large coefficients
def largecoef(rng):
    for k in range(25):
        digits = rng.choice([10, 20, 40, 60])
        lo, hi = -10**digits, 10**digits
        nf = rng.randint(1, 3)
        den = 1
        for _ in range(nf):
            dg = rng.randint(1, 3)
            while True:
                p = rpoly(rng, dg, lo, hi)
                if dg == 1 or irreducible(p):
                    break
            den *= p ** rng.randint(1, 2)
        dd = Poly(den, x).degree()
        f = expand(rpoly(rng, max(dd - 1, 0), lo, hi)) / expand(den)
        add("largecoef", f, f"{digits}-digit coefficients, deg {dd}")


# ---------------------------------------------------------------- stress / blowup
def stress(rng):
    for n in range(4, 11):
        add("stress", 1/prod([x**2 + k for k in range(1, n + 1)]), f"prod_(k<={n}) (x^2+k)")
    sd3 = expand(prod([x + s2*sqrt(2) + s3*sqrt(3) + s5*sqrt(5)
                       for s2 in (1, -1) for s3 in (1, -1) for s5 in (1, -1)]))
    add("stress", 1/sd3, "Swinnerton-Dyer deg 8")
    add("stress", x**3/sd3, "Swinnerton-Dyer deg 8")
    for n in (20, 24, 30):
        add("stress", 1/(x**n + 1), f"1/(x^{n}+1)")
    add("stress", 1/(x**24 - 1), "1/(x^24-1)")
    for k in (2, 3):
        add("stress", 1/(x**7 + x + 1)**k, "irreducible deg 7, multiplicity")
    add("stress", 1/((x**3 + x + 1)**3*(x**2 + 2)**2), "mixed multiplicities")
    for dg in (12, 16, 20):
        d = rirred(rng, dg, -2, 2)
        add("stress", rnum(rng, dg - 1, -2, 2)/d, f"irreducible deg {dg}")
    add("stress", 1/(x**10 - 2*(10*x - 1)**2), "Mignotte-like")
    add("stress", 1/(x**8 + x**4 + 1)**3, "Hermite deg 24")
    add("stress", (x**6 + 1)/(x**12 + x**6 + 1), "palindromic deg 12")
    add("stress", 1/((x**2 + 1)**5 * (x**2 + 2)**5), "Hermite two quadratics")
    add("stress", 1/(x**9 + 3*x**3 + 1), "")
    add("stress", 1/(x**6 - 2*x**3 - 1)**2, "")


# ---------------------------------------------------------------- Czichowski-favourable
def czichowski(rng):
    """High-degree squarefree denominators whose log part has few distinct
    residues: the subresultant chain of LRT is long, the resultant has small
    factors, and the lex Groebner basis of Czichowski is tiny (Bronstein §2.6).
    Fixed seed of its own, so the earlier families keep their ids."""
    rng = random.Random(2661)
    k = 0
    while k < 8:
        S1 = rirred(rng, rng.randint(4, 7), -3, 3)
        S2 = rirred(rng, rng.randint(3, 6), -3, 3)
        if expand(S1 - S2) == 0:
            continue
        c1, c2 = rng.choice([1, 2, 3, -1, -2]), rng.choice([1, -1, 3, -3, Rational(1, 2)])
        f = cancel(together(c1*diff(S1, x)/S1 + c2*diff(S2, x)/S2))
        if Poly(f.as_numer_denom()[1], x).degree() < 7:
            continue
        add("czichowski", f, f"{c1} log S1 + {c2} log S2, deg {Poly(S1, x).degree()}+{Poly(S2, x).degree()}")
        k += 1
    k = 0
    while k < 6:
        P = rpoly(rng, rng.randint(5, 7), -3, 3)
        Q = rpoly(rng, rng.randint(3, 6), -3, 3)
        f = cancel(together((diff(P, x)*Q - P*diff(Q, x)) / (P**2 + Q**2)))
        if f == 0:
            continue
        add("czichowski", f, f"D[ArcTan[P/Q]], deg P = {Poly(P, x).degree()}, deg Q = {Poly(Q, x).degree()}")
        k += 1


# ---------------------------------------------------------------- LogToReal-relevant
def logtoreal(rng):
    """Resultant factors of degree 3-4 with complex roots, beyond the shape
    dispatch: quartics irreducible over Q that split over Q(Sqrt[d]) into two
    quadratics with complex roots (Re and Im of the residues are quadratic
    surds: general LogToReal in radicals), and random cubics and quartics
    (Re and Im of higher degree: RootSum, or Root objects)."""
    rng = random.Random(2662)
    k = 0
    while k < 10:
        a_, b_, c_, e_ = rng.randint(-3, 3), rng.randint(2, 6), rng.randint(1, 2), rng.randint(-2, 2)
        d_ = rng.choice([2, 3, 5, 6, 7])
        D = expand((x**2 + a_*x + b_)**2 - d_*(c_*x + e_)**2)
        if not irreducible(D):
            continue
        roots_real = [r for r in Poly(D, x).nroots() if abs(r.as_real_imag()[1]) < 1e-12]
        if roots_real:
            continue
        add("logtoreal", rnum(rng, 3, -4, 4) / D, f"quartic split over Q(Sqrt[{d_}])")
        k += 1
    k = 0
    while k < 5:
        D = x**3 + rng.randint(-4, 4)*x**2 + rng.randint(-4, 4)*x + rng.choice([1, 2, 3, -1, -2, 5])
        if not irreducible(D):
            continue
        add("logtoreal", rnum(rng, 2, -3, 3) / D, "irreducible cubic")
        k += 1
    k = 0
    while k < 5:
        D = rirred(rng, 4, -3, 3)
        add("logtoreal", rnum(rng, 3, -3, 3) / D, "irreducible quartic")
        k += 1


def main():
    rng = random.Random(266)
    classical(); bronstein(); sqfree(rng); nonsqfree(rng); multiplicity(); rootsum(rng)
    rioboo(rng); param_generic(); param_special(); alg_qi(rng); alg_sqrt(rng); alg_nested()
    alg_root(); largecoef(rng); stress(rng)
    czichowski(rng); logtoreal(rng)
    out = os.path.join(HERE, "corpus.tsv")
    with open(out, "w") as fh:
        fh.write("# id\tfamily\tintegrand\tsympy\tnote\n")
        for row in CASES:
            fh.write("\t".join(row) + "\n")
    fams = {}
    for row in CASES:
        fams[row[1]] = fams.get(row[1], 0) + 1
    print(f"wrote {out}: {len(CASES)} cases", fams, file=sys.stderr)


if __name__ == "__main__":
    main()
