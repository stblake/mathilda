#!/usr/bin/env python3
"""Experiment 06 -- Multivariate factorisation and expansion (sympy column).

Same six kernels as ``factor_multivariate.m``, same order.  Checks are term
counts and factor counts, both ordering-independent.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sympy

from harness import bench, check, require
from data import multi_poly

require(["sympy", "sympy:factor", "sympy:factor_list", "sympy:expand"])

x, y, z, u, v = sympy.symbols("x y z u v")


def nterms(e):
    """Length[] of an expanded Plus == its term count."""
    return len(sympy.Add.make_args(sympy.expand(e)))


def nfactors(e):
    coeff, factors = sympy.factor_list(e)
    return len(factors) + 1


bench("Expand (x+y+z+i), 6 factors", lambda: multi_poly(6, [x, y, z]))
check("Expand (x+y+z+i), 6 factors", nterms(multi_poly(6, [x, y, z])))

m6 = multi_poly(6, [x, y, z])
bench("Factor 3-var, 6 factors", lambda: sympy.factor(m6))
check("Factor 3-var, 6 factors", nfactors(m6))

m5 = multi_poly(4, [x, y, z, u, v])
bench("Factor 5-var, 4 factors", lambda: sympy.factor(m5))
check("Factor 5-var, 4 factors", nfactors(m5))

irr = sympy.expand(x ** 3 + y ** 3 + z ** 3 - 3 * x * y * z + 1)
bench("Factor irreducible 3-var cubic", lambda: sympy.factor(irr))
check("Factor irreducible 3-var cubic", nfactors(irr))

bench("Expand (x+y+z)^12", lambda: sympy.expand((x + y + z) ** 12))
check("Expand (x+y+z)^12", nterms((x + y + z) ** 12))

d2 = x ** 24 - y ** 24
bench("Factor x^24 - y^24", lambda: sympy.factor(d2))
check("Factor x^24 - y^24", nfactors(d2))
