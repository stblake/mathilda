#!/usr/bin/env python3
"""Experiment 14 -- Numerical root finding (scipy.optimize column).

Same five kernels as ``root_finding.m``, same order.  Checks are the root value
rounded, or the root count -- both invariant under the iteration path taken.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import numpy as np
from scipy.optimize import brentq, root

from harness import bench, check, require
from data import dense_upoly

require(["scipy.optimize:brentq", "scipy.optimize:root", "numpy:roots"])


def r(x, p):
    return int(math.floor(x * 10 ** p + 0.5))


f1 = lambda t: math.cos(t) - t
bench("FindRoot Cos[x]==x", lambda: brentq(f1, 0.0, 2.0, xtol=1e-14))
check("FindRoot Cos[x]==x", r(brentq(f1, 0.0, 2.0, xtol=1e-14), 6))

f2 = lambda t: t ** 3 - 2 * t - 5
bench("FindRoot x^3-2x-5", lambda: brentq(f2, 1.0, 3.0, xtol=1e-14))
check("FindRoot x^3-2x-5", r(brentq(f2, 1.0, 3.0, xtol=1e-14), 6))

f3 = lambda v: [v[0] ** 2 + v[1] ** 2 - 4.0, v[0] * v[1] - 1.0]
bench("FindRoot 2x2 nonlinear system",
      lambda: root(f3, [1.8, 0.6], tol=1e-13))
check("FindRoot 2x2 nonlinear system",
      r(root(f3, [1.8, 0.6], tol=1e-13).x[0], 6))

# numpy.roots takes coefficients highest-degree first.
import sympy
_x = sympy.Symbol("x")
_coeffs = [float(c) for c in sympy.Poly(dense_upoly(50, _x), _x).all_coeffs()]
bench("NSolve all roots, degree 50", lambda: np.roots(_coeffs))
check("NSolve all roots, degree 50", len(np.roots(_coeffs)))

f5 = lambda t: math.tan(t) - t
# Bracket inside one branch of tan, away from the pole at pi/2 * 3.
bench("FindRoot Tan[x]==x near 4.5",
      lambda: brentq(f5, 4.0, 4.7, xtol=1e-13))
check("FindRoot Tan[x]==x near 4.5", r(brentq(f5, 4.0, 4.7, xtol=1e-13), 4))
