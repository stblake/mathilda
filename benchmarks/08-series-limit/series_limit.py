#!/usr/bin/env python3
"""Experiment 08 -- Series expansion and limits (sympy column).

Same six kernels as ``series_limit.m``, same order.

Series and Limit fail in opposite directions, which is why both are here.
Series is a cost problem (the answer is never in doubt).  Limit is a correctness
problem: the hard cases need a comparability ordering, and a system without one
returns the input unevaluated -- fast and wrong.  Checks are single coefficients
and exact limit values.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sympy

from harness import bench, check, require

require(["sympy", "sympy:series", "sympy:limit", "sympy:oo"])

x = sympy.Symbol("x")


def coeff(expr, n, order):
    """The x**n coefficient of a series to `order` -- matches SeriesCoefficient."""
    return sympy.series(expr, x, 0, order + 1).removeO().coeff(x, n)


bench("Series Exp[Sin[x]] to order 20",
      lambda: sympy.series(sympy.exp(sympy.sin(x)), x, 0, 21))
check("Series Exp[Sin[x]] to order 20",
      coeff(sympy.exp(sympy.sin(x)), 12, 20))

bench("Series 1/(1-x-x^2) to order 60",
      lambda: sympy.series(1 / (1 - x - x ** 2), x, 0, 61))
check("Series 1/(1-x-x^2) to order 60",
      coeff(1 / (1 - x - x ** 2), 50, 60))

bench("Series Log[1+Sin[x]] to order 24",
      lambda: sympy.series(sympy.log(1 + sympy.sin(x)), x, 0, 25))
check("Series Log[1+Sin[x]] to order 24",
      coeff(sympy.log(1 + sympy.sin(x)), 15, 24))

bench("Limit Sin[x]/x at 0", lambda: sympy.limit(sympy.sin(x) / x, x, 0))
check("Limit Sin[x]/x at 0", sympy.limit(sympy.sin(x) / x, x, 0))

e5 = (sympy.exp(x) - 1 - x) / x ** 2
bench("Limit (Exp[x]-1-x)/x^2 at 0", lambda: sympy.limit(e5, x, 0))
check("Limit (Exp[x]-1-x)/x^2 at 0", sympy.limit(e5, x, 0))

e6 = x * (sympy.log(x + 1) - sympy.log(x))
bench("Limit x (Log[x+1]-Log[x]) at Infinity",
      lambda: sympy.limit(e6, x, sympy.oo))
check("Limit x (Log[x+1]-Log[x]) at Infinity", sympy.limit(e6, x, sympy.oo))
