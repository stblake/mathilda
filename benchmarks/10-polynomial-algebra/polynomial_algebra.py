#!/usr/bin/env python3
"""Experiment 10 -- Polynomial algebra: GCD, resultants, canonical forms (sympy).

Same six kernels as ``polynomial_algebra.m``, same order.

WHY IT IS A SEPARATE EXPERIMENT.  Integrate, Simplify, Solve and Factor all call
these primitives.  When one of those rows is slow, this experiment answers
whether the cost is in the high-level algorithm or in the arithmetic underneath.
Checks are degrees and evaluated values, both ordering-independent.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sympy

from harness import bench, check, require
from data import dense_upoly, factorable_upoly

require(["sympy", "sympy:gcd", "sympy:resultant", "sympy:discriminant",
         "sympy:cancel", "sympy:div"])

x = sympy.Symbol("x")

sh = factorable_upoly(10, x)
a1 = sympy.expand(sh * dense_upoly(40, x))
b1 = sympy.expand(sh * dense_upoly(38, x))

bench("PolynomialGCD, shared deg-20 factor", lambda: sympy.gcd(a1, b1))
check("PolynomialGCD, shared deg-20 factor", sympy.degree(sympy.gcd(a1, b1), x))

c1, c2 = dense_upoly(40, x), dense_upoly(39, x) + 1
bench("PolynomialGCD, coprime deg 40", lambda: sympy.gcd(c1, c2))
check("PolynomialGCD, coprime deg 40", sympy.degree(sympy.gcd(c1, c2), x))

d20 = dense_upoly(20, x)
bench("Discriminant of deg 20", lambda: sympy.discriminant(d20, x))
check("Discriminant of deg 20", int(sympy.discriminant(d20, x)) % 1000003)

bench("Cancel deg-60 over deg-58", lambda: sympy.cancel(a1 / b1))
check("Cancel deg-60 over deg-58",
      int(sympy.floor(sympy.N(sympy.cancel(a1 / b1).subs(
          x, sympy.Rational(11, 10)) * 10 ** 6, 50) + sympy.Rational(1, 2))))

bench("PolynomialQuotient deg 60 / deg 20",
      lambda: sympy.div(a1, sympy.expand(sh), x)[0])
check("PolynomialQuotient deg 60 / deg 20",
      sympy.degree(sympy.div(a1, sympy.expand(sh), x)[0], x))

bench("Expand (1+x)^400", lambda: sympy.expand((1 + x) ** 400))
check("Expand (1+x)^400", sympy.expand((1 + x) ** 400).coeff(x, 200))
