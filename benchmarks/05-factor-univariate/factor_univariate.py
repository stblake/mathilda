#!/usr/bin/env python3
"""Experiment 05 -- Univariate polynomial factorisation (sympy column).

Same six kernels as ``factor_univariate.m``, same order.

The Mathilda side routes through FLINT when USE_FLINT is compiled in; sympy uses
its own pure-Python Zassenhaus.  Check values are FACTOR COUNTS: output ordering
is not canonical across systems, but the number of irreducible factors is.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sympy

from harness import bench, check, require
from data import dense_upoly, factorable_upoly, sparse_upoly

require(["sympy", "sympy:factor", "sympy:factor_list", "sympy:gcd"])

x = sympy.Symbol("x")


def nfactors(e):
    """Match Length[FactorList[...]]: the unit term plus each distinct factor."""
    coeff, factors = sympy.factor_list(e)
    return len(factors) + 1


f1 = sympy.expand(factorable_upoly(8, x))
bench("Factor product of 8 quadratics", lambda: sympy.factor(f1))
check("Factor product of 8 quadratics", nfactors(f1))

f2 = sympy.expand(factorable_upoly(16, x))
bench("Factor product of 16 quadratics", lambda: sympy.factor(f2))
check("Factor product of 16 quadratics", nfactors(f2))

f3 = x ** 120 - 1
bench("Factor x^120 - 1", lambda: sympy.factor(f3))
check("Factor x^120 - 1", nfactors(f3))

f4 = sparse_upoly(60, x)
bench("Factor sparse degree 60", lambda: sympy.factor(f4))
check("Factor sparse degree 60", nfactors(f4))

f5 = dense_upoly(60, x)
bench("Factor dense degree 60", lambda: sympy.factor(f5))
check("Factor dense degree 60", nfactors(f5))

g1 = sympy.expand(factorable_upoly(10, x))
g2 = sympy.expand(factorable_upoly(6, x) * (x ** 4 + x + 1))
bench("PolynomialGCD of two deg-20 polys", lambda: sympy.gcd(g1, g2))
check("PolynomialGCD of two deg-20 polys",
      sympy.degree(sympy.gcd(g1, g2), x))
