#!/usr/bin/env python3
"""Experiment 125 -- definite integrals by Meijer G functions (sympy/mpmath column).

Same cases, labels and checks as ``meijer_g_integrals.m``; see its header.

* ``roach ...``  sympy ``integrate((f, (x, 0, oo)), meijerg=True)`` -- meijerint,
  the same Meijer G convolution method (PBM 2.24.1) as Mathilda's
  ``Method -> "Roach"``.  sympy caches every intermediate, so each timed call
  starts from a cleared cache; otherwise runs 2-4 would time a dictionary
  lookup.
* ``hyperexpand ...``  ``sympy.hyperexpand`` against ``FunctionExpand`` of the
  same pFq (Roach's ISSAC '96 gallery).
* ``N ...``  30-digit numerics, ``mpmath.hyper`` / ``mpmath.meijerg``.

A case sympy leaves unevaluated checks as ``unevaluated``: numericising a
returned ``Integral`` would quietly integrate it numerically and pass the gate.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sympy
import mpmath
from sympy import (Symbol, Rational as R, oo, exp, sqrt, sin, cos, log, pi,
                   besselj, besseli, erfc, Si, hyper, hyperexpand, Integral)
from sympy.core.cache import clear_cache

from harness import bench, check, require

require(["sympy", "sympy:integrate", "sympy:hyperexpand",
         "sympy.integrals.meijerint:meijerint_definite",
         "mpmath:hyper", "mpmath:meijerg"])

x = Symbol("x", positive=True)
z = Symbol("z")


def meij(f):
    clear_cache()
    return sympy.integrate(f, (x, 0, oo), meijerg=True)


def chk(v, at=None):
    """Round[10^6 value], or 'unevaluated' when sympy kept an integral/pFq."""
    if v.has(Integral) or v.has(hyper) or v.has(sympy.meijerg):
        return '"unevaluated"'   # the .m side prints InputForm of a String
    if at is not None:
        v = v.subs(z, at)
    return int(sympy.floor(sympy.re(sympy.N(v * 10**6, 40)) + R(1, 2)))


INTEGRALS = [
    ("roach Exp[-3x/2] BesselJ[0,x]^2", exp(-R(3, 2) * x) * besselj(0, x)**2),
    ("roach x^(-1/4)/Sqrt[x^4+1]", x**R(-1, 4) / sqrt(x**4 + 1)),
    ("roach x^(-1/2) (1+2x)^(-3/2)", x**R(-1, 2) * (1 + 2 * x)**R(-3, 2)),
    ("roach x^(-1/2) (3+x)^-2", x**R(-1, 2) * (3 + x)**-2),
    ("roach BesselJ[1,x] Exp[-x^2/2]", besselj(1, x) * exp(-x**2 / 2)),
    ("roach Exp[-x] Sqrt[x]/(2+x)", exp(-x) * sqrt(x) / (2 + x)),
    ("roach x Erfc[x] Exp[-x^2/2]", x * erfc(x) * exp(-x**2 / 2)),
    ("roach (Si[2x]-Pi/2) Exp[-x]", (Si(2 * x) - pi / 2) * exp(-x)),
    ("roach Sin[x]/(x^2+4)^(3/4)", sin(x) / (x**2 + 4)**R(3, 4)),
    ("roach BesselJ[1,x] BesselJ[1,2x]/x^2", besselj(1, x) * besselj(1, 2 * x) / x**2),
    ("roach x^2 Exp[-x^2] BesselJ[0,x]", x**2 * exp(-x**2) * besselj(0, x)),
    ("roach x^(-1/2) Exp[-x] BesselJ[0,x]", x**R(-1, 2) * exp(-x) * besselj(0, x)),
]

for label, f in INTEGRALS:
    bench(label, lambda f=f: meij(f))
    check(label, chk(meij(f)))

PFQS = [
    ("hyperexpand 1F1(1/2; 9/2; z)", hyper([R(1, 2)], [R(9, 2)], z)),
    ("hyperexpand 2F1(-3/2,-1/2; 2; z)", hyper([R(-3, 2), R(-1, 2)], [2], z)),
    ("hyperexpand 3F2(1,2,3; 1/2,4; z)", hyper([1, 2, 3], [R(1, 2), 4], z)),
    ("hyperexpand 1F2(3/2; 5/2,5; z)", hyper([R(3, 2)], [R(5, 2), 5], z)),
]


def hexp(h):
    clear_cache()
    return hyperexpand(h)


for label, h in PFQS:
    bench(label, lambda h=h: hexp(h))
    check(label, chk(hexp(h), R(3, 10)))

# ---- 30-digit numerics -------------------------------------------------------
mpmath.mp.dps = 30


def n_pfq():
    return mpmath.hyper([mpmath.mpf(-15) / 4, 1, 1], [mpmath.mpf(-7) / 2, 2], 1)


def n_meijer():
    # G^{1,4}_{4,2}(4 | 1/2,1,1,1 ; 1, 0), the Roach form of Log[1+x] Cos[x]/x
    return mpmath.meijerg([[0.5, 1, 1, 1], []], [[1], [0]], 4)


bench("N 3F2(-15/4,1,1; -7/2,2; 1) 30 digits", n_pfq)
check("N 3F2(-15/4,1,1; -7/2,2; 1) 30 digits", int(mpmath.floor(n_pfq() * 10**6 + 0.5)))
bench("N MeijerG G^{1,4}_{4,2}(4) 30 digits", n_meijer)
check("N MeijerG G^{1,4}_{4,2}(4) 30 digits", int(mpmath.floor(n_meijer() * 10**6 + 0.5)))
