#!/usr/bin/env python3
"""Experiment 02 -- Symbolic integration, transcendental cases (sympy column).

Same six kernels as ``integrate_transcendental.m``, same order.

Each row exercises a different branch of the Risch machinery, so a slow row
names a subsystem rather than "integration": exp tower, log tower, mixed
(coupled) tower, tan tower, a rational function *of* an exponential, and one
genuinely non-elementary integrand.

THE LAST ROW IS THE POINT.  ``exp(x)/x`` has no elementary antiderivative.  A
correct system proves that and returns Ei; an incorrect one returns the integral
unevaluated -- which is FAST and WRONG, and would read as a win without a value
check.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sympy

from harness import bench, check, require

require(["sympy", "sympy:integrate", "sympy:Ei"])

x = sympy.Symbol("x")


def dif(expr, a, b):
    """Round[10**6 (F(b) - F(a))] -- invariant under +C and under regrouping."""
    F = sympy.integrate(expr, x)
    d = F.subs(x, b) - F.subs(x, a)
    return int(sympy.floor(sympy.re(sympy.N(d * 10 ** 6, 50)) + sympy.Rational(1, 2)))


R = sympy.Rational

bench("integrate x Exp[x^2]", lambda: sympy.integrate(x * sympy.exp(x ** 2), x))
check("integrate x Exp[x^2]", dif(x * sympy.exp(x ** 2), R(1, 4), R(3, 4)))

bench("integrate Log[x]^3", lambda: sympy.integrate(sympy.log(x) ** 3, x))
check("integrate Log[x]^3", dif(sympy.log(x) ** 3, R(3, 2), R(5, 2)))

bench("integrate Sin[x] Exp[x]",
      lambda: sympy.integrate(sympy.sin(x) * sympy.exp(x), x))
check("integrate Sin[x] Exp[x]",
      dif(sympy.sin(x) * sympy.exp(x), R(1, 2), R(3, 2)))

bench("integrate Tan[x]^3", lambda: sympy.integrate(sympy.tan(x) ** 3, x))
check("integrate Tan[x]^3", dif(sympy.tan(x) ** 3, R(1, 4), R(1, 2)))

bench("integrate 1/(1+Exp[x])",
      lambda: sympy.integrate(1 / (1 + sympy.exp(x)), x))
check("integrate 1/(1+Exp[x])", dif(1 / (1 + sympy.exp(x)), R(1, 2), R(3, 2)))

# Non-elementary: must return Ei rather than give up.
bench("integrate Exp[x]/x (non-elementary)",
      lambda: sympy.integrate(sympy.exp(x) / x, x))
check("integrate Exp[x]/x (non-elementary)",
      dif(sympy.exp(x) / x, R(3, 2), R(5, 2)))
