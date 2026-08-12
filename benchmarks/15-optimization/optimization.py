#!/usr/bin/env python3
"""Experiment 15 -- Numerical optimisation (scipy.optimize column).

Same five kernels as ``optimization.m``, same order.  Checks are the objective
value at the minimum, rounded -- invariant under which descent path was taken.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import numpy as np
from scipy.optimize import minimize

from harness import bench, check, require

require(["scipy.optimize:minimize", "numpy.polynomial:polynomial",
         "numpy:polyfit"])

TOL = dict(tol=1e-12)


def r(x, p):
    return int(math.floor(x * 10 ** p + 0.5))


q = lambda v: (v[0] - 3.0) ** 2 + 2.0
bench("FindMinimum quadratic", lambda: minimize(q, [0.0], **TOL))
check("FindMinimum quadratic", r(minimize(q, [0.0], **TOL).fun, 6))

ros = lambda v: (1.0 - v[0]) ** 2 + 100.0 * (v[1] - v[0] ** 2) ** 2
bench("FindMinimum Rosenbrock 2-D", lambda: minimize(ros, [-1.0, 1.0], **TOL))
check("FindMinimum Rosenbrock 2-D", r(minimize(ros, [-1.0, 1.0], **TOL).fun, 4))

xs = lambda v: v[0] * math.sin(v[0])
bench("FindMinimum x Sin[x] on [2,6]", lambda: minimize(xs, [4.0], **TOL))
check("FindMinimum x Sin[x] on [2,6]", r(minimize(xs, [4.0], **TOL).fun, 4))

ti = np.arange(1.0, 201.0)
yi = 2.0 * ti + 1.0
bench("Fit linear, 200 points", lambda: np.polyfit(ti, yi, 1))
_c1 = np.polyfit(ti, yi, 1)
check("Fit linear, 200 points", r(float(np.polyval(_c1, 10.0)), 6))

t2 = np.arange(1, 501) / 100.0
y2 = t2 ** 3 - 2.0 * t2 + 1.0
bench("Fit degree-5 polynomial, 500 points", lambda: np.polyfit(t2, y2, 5))
_c2 = np.polyfit(t2, y2, 5)
check("Fit degree-5 polynomial, 500 points", r(float(np.polyval(_c2, 2.0)), 4))
