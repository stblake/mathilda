#!/usr/bin/env python3
"""Experiment 65 -- L-BFGS-B accuracy + efficiency stress (scipy column).

Same 10 hard test functions as ``lbfgsb_accuracy.m``, same labels, same order,
each with its analytic Jacobian passed to scipy.optimize.minimize(method=
"L-BFGS-B") so the comparison is solver-vs-solver. The check is the objective
at the optimum rounded to 6 digits (0 for the f*=0 problems) or the exact
integer optimum (Trid -50, bound corner 5): reaching it is the accuracy gate.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import minimize

from harness import bench, check, require

require(["scipy.optimize:minimize"])

_OPTS = {"maxiter": 15000, "maxfun": 200000, "ftol": 1e-15, "gtol": 1e-10}


def _solve(fun, jac, x0, bounds=None):
    return minimize(fun, np.asarray(x0, float), jac=jac, method="L-BFGS-B",
                    bounds=bounds, options=_OPTS)


# ---- 01 Rosenbrock 2D, f* = 0 ----
def _rb(v):   x, y = v; return (1 - x) ** 2 + 100 * (y - x * x) ** 2
def _rbg(v):  x, y = v; return np.array([-2 * (1 - x) - 400 * x * (y - x * x), 200 * (y - x * x)])

# ---- ill-conditioned diagonal quadratic ----
def _ill(n, c):
    w = 10.0 ** (c * np.arange(n) / (n - 1))
    return (lambda v: float(np.sum(w * v * v)), lambda v: 2.0 * w * v, np.ones(n))

# ---- 03 Beale, f* = 0 at (3, 0.5) ----
def _beale(v):
    x, y = v
    return (1.5 - x + x * y) ** 2 + (2.25 - x + x * y * y) ** 2 + (2.625 - x + x * y ** 3) ** 2
def _bealeg(v):
    x, y = v
    a, b, c = 1.5 - x + x * y, 2.25 - x + x * y * y, 2.625 - x + x * y ** 3
    return np.array([2 * a * (y - 1) + 2 * b * (y * y - 1) + 2 * c * (y ** 3 - 1),
                     2 * a * x + 4 * b * x * y + 6 * c * x * y * y])

# ---- 04 Booth, f* = 0 at (1, 3) ----
def _booth(v):  x, y = v; return (x + 2 * y - 7) ** 2 + (2 * x + y - 5) ** 2
def _boothg(v):
    x, y = v
    return np.array([2 * (x + 2 * y - 7) + 4 * (2 * x + y - 5),
                     4 * (x + 2 * y - 7) + 2 * (2 * x + y - 5)])

# ---- 05 Matyas, f* = 0 at origin ----
def _mat(v):   x, y = v; return 0.26 * (x * x + y * y) - 0.48 * x * y
def _matg(v):  x, y = v; return np.array([0.52 * x - 0.48 * y, 0.52 * y - 0.48 * x])

# ---- 06 Wood 4D, f* = 0 (order x, y, z, w) ----
def _wood(v):
    x, y, z, w = v
    return (100 * (y - x * x) ** 2 + (1 - x) ** 2 + 90 * (w - z * z) ** 2 + (1 - z) ** 2
            + 10.1 * ((y - 1) ** 2 + (w - 1) ** 2) + 19.8 * (y - 1) * (w - 1))
def _woodg(v):
    x, y, z, w = v
    return np.array([-400 * x * (y - x * x) - 2 * (1 - x),
                     200 * (y - x * x) + 20.2 * (y - 1) + 19.8 * (w - 1),
                     -360 * z * (w - z * z) - 2 * (1 - z),
                     180 * (w - z * z) + 20.2 * (w - 1) + 19.8 * (y - 1)])

# ---- 07 Powell singular 4D, f* = 0 (order x, y, z, w) ----
def _pow(v):
    x, y, z, w = v
    return (x + 10 * y) ** 2 + 5 * (z - w) ** 2 + (y - 2 * z) ** 4 + 10 * (x - w) ** 4
def _powg(v):
    x, y, z, w = v
    return np.array([2 * (x + 10 * y) + 40 * (x - w) ** 3,
                     20 * (x + 10 * y) + 4 * (y - 2 * z) ** 3,
                     10 * (z - w) - 8 * (y - 2 * z) ** 3,
                     -10 * (z - w) - 40 * (x - w) ** 3])

# ---- 08 Trid n, f* = -n(n+4)(n-1)/6 ----
def _trid(v):
    return float(np.sum((v - 1) ** 2) - np.sum(v[1:] * v[:-1]))
def _tridg(v):
    n = v.size
    g = 2 * (v - 1)
    g[1:] -= v[:-1]
    g[:-1] -= v[1:]
    return g

# ---- 09 bound-active quadratic on [0,1]^2, f* = 5 ----
def _bq(v):   x, y = v; return (x - 2) ** 2 + (y - 3) ** 2
def _bqg(v):  x, y = v; return np.array([2 * (x - 2), 2 * (y - 3)])


ill10f, ill10g, ill10x = _ill(10, 6)
ill50f, ill50g, ill50x = _ill(50, 4)

bench("01 Rosenbrock2D", lambda: _solve(_rb, _rbg, [-1.2, 1.0]))
check("01 Rosenbrock2D", int(np.floor(1e6 * _solve(_rb, _rbg, [-1.2, 1.0]).fun + 0.5)))

bench("02 illcond n10 c1e6", lambda: _solve(ill10f, ill10g, ill10x))
check("02 illcond n10 c1e6", int(np.floor(1e6 * _solve(ill10f, ill10g, ill10x).fun + 0.5)))

bench("03 Beale", lambda: _solve(_beale, _bealeg, [1.0, 1.0]))
check("03 Beale", int(np.floor(1e6 * _solve(_beale, _bealeg, [1.0, 1.0]).fun + 0.5)))

bench("04 Booth", lambda: _solve(_booth, _boothg, [0.0, 0.0]))
check("04 Booth", int(np.floor(1e6 * _solve(_booth, _boothg, [0.0, 0.0]).fun + 0.5)))

bench("05 Matyas", lambda: _solve(_mat, _matg, [3.0, 3.0]))
check("05 Matyas", int(np.floor(1e6 * _solve(_mat, _matg, [3.0, 3.0]).fun + 0.5)))

bench("06 Wood", lambda: _solve(_wood, _woodg, [-3.0, -1.0, -3.0, -1.0]))
check("06 Wood", int(np.floor(1e6 * _solve(_wood, _woodg, [-3.0, -1.0, -3.0, -1.0]).fun + 0.5)))

bench("07 Powell-singular", lambda: _solve(_pow, _powg, [3.0, -1.0, 0.0, 1.0]))
check("07 Powell-singular", int(np.floor(1e6 * _solve(_pow, _powg, [3.0, -1.0, 0.0, 1.0]).fun + 0.5)))

bench("08 Trid n6", lambda: _solve(_trid, _tridg, np.zeros(6)))
check("08 Trid n6", int(np.floor(_solve(_trid, _tridg, np.zeros(6)).fun + 0.5)))

bench("09 Bound-corner", lambda: _solve(_bq, _bqg, [0.0, 0.0], bounds=[(0.0, 1.0), (0.0, 1.0)]))
check("09 Bound-corner", int(np.floor(_solve(_bq, _bqg, [0.0, 0.0], bounds=[(0.0, 1.0), (0.0, 1.0)]).fun + 0.5)))

bench("10 illcond n50 c1e4", lambda: _solve(ill50f, ill50g, ill50x))
check("10 illcond n50 c1e4", int(np.floor(1e6 * _solve(ill50f, ill50g, ill50x).fun + 0.5)))
