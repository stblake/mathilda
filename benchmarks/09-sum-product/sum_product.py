#!/usr/bin/env python3
"""Experiment 09 -- Symbolic summation and products (sympy column).

Same six kernels as ``sum_product.m``, same order.

Symbolic summation is where a CAS visibly either has an algorithm or does not:
Gosper either finds the hypergeometric antidifference or proves none exists.
There is no slow-but-correct middle, so these rows split cleanly into "fast" and
"absent".  Checks are exact closed-form values at a small n.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sympy

from harness import bench, check, require

require(["sympy", "sympy:Sum", "sympy:Product", "sympy:binomial", "sympy:zeta"])

k, n, r = sympy.symbols("k n r")


def dosum(term, lo, hi):
    return sympy.Sum(term, (k, lo, hi)).doit()


def doprod(term, lo, hi):
    return sympy.Product(term, (k, lo, hi)).doit()


bench("Sum k^5 to n, closed form", lambda: dosum(k ** 5, 1, n))
check("Sum k^5 to n, closed form",
      sympy.simplify(dosum(k ** 5, 1, n).subs(n, 10)))

bench("Sum 1/k^2 to Infinity", lambda: dosum(1 / k ** 2, 1, sympy.oo))
check("Sum 1/k^2 to Infinity",
      int(sympy.floor(sympy.N(dosum(1 / k ** 2, 1, sympy.oo) * 10 ** 6, 50)
                      + sympy.Rational(1, 2))))

bench("Sum r^k to n, symbolic ratio", lambda: dosum(r ** k, 0, n))
check("Sum r^k to n, symbolic ratio",
      sympy.simplify(dosum(r ** k, 0, n).subs({r: 3, n: 6})))

bench("Sum Binomial[n,k] over k", lambda: dosum(sympy.binomial(n, k), 0, n))
check("Sum Binomial[n,k] over k",
      sympy.simplify(dosum(sympy.binomial(n, k), 0, n).subs(n, 12)))

bench("Sum 1/(k(k+1)) to n", lambda: dosum(1 / (k * (k + 1)), 1, n))
check("Sum 1/(k(k+1)) to n",
      int(sympy.floor(sympy.N(dosum(1 / (k * (k + 1)), 1, n).subs(n, 10)
                              * 10 ** 6, 50) + sympy.Rational(1, 2))))

bench("Product (1+1/k) to n", lambda: doprod(1 + 1 / k, 1, n))
check("Product (1+1/k) to n",
      sympy.simplify(doprod(1 + 1 / k, 1, n).subs(n, 10)))
