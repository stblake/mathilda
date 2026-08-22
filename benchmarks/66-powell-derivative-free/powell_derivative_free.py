#!/usr/bin/env python3
"""Experiment 66 -- Powell derivative-free optimization (scipy column).

Same 13 test functions as ``powell_derivative_free.m``, same labels, same order,
solved with scipy.optimize.minimize(method="Powell"). Unlike experiment 65
(L-BFGS-B), NO Jacobian is passed -- Powell is derivative-free, so the fair race
is derivative-free-vs-derivative-free. The check is the objective at the optimum
rounded to 6 digits (0 for the f*=0 problems) or the exact integer optimum
(Trid -50, bound corner 5): reaching it is the accuracy gate.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import minimize

from harness import bench, check, require

require(["scipy.optimize:minimize"])

_OPTS = {"maxiter": 20000, "maxfev": 200000, "xtol": 1e-10, "ftol": 1e-12}


def _solve(fun, x0, bounds=None):
    return minimize(fun, np.asarray(x0, float), method="Powell",
                    bounds=bounds, options=_OPTS)


# ---- 01 Rosenbrock 2D, f* = 0 ----
def _rb(v):   x, y = v; return (1 - x) ** 2 + 100 * (y - x * x) ** 2

# ---- 02 Beale, f* = 0 at (3, 0.5) ----
def _beale(v):
    x, y = v
    return (1.5 - x + x * y) ** 2 + (2.25 - x + x * y * y) ** 2 + (2.625 - x + x * y ** 3) ** 2

# ---- 03 Booth, f* = 0 at (1, 3) ----
def _booth(v):  x, y = v; return (x + 2 * y - 7) ** 2 + (2 * x + y - 5) ** 2

# ---- 04 Matyas, f* = 0 at origin ----
def _mat(v):   x, y = v; return 0.26 * (x * x + y * y) - 0.48 * x * y

# ---- 05 Powell singular 4D, f* = 0 (order x, y, z, w) ----
def _pow(v):
    x, y, z, w = v
    return (x + 10 * y) ** 2 + 5 * (z - w) ** 2 + (y - 2 * z) ** 4 + 10 * (x - w) ** 4

# ---- 06 Wood 4D, f* = 0 (order x, y, z, w) ----
def _wood(v):
    x, y, z, w = v
    return (100 * (y - x * x) ** 2 + (1 - x) ** 2 + 90 * (w - z * z) ** 2 + (1 - z) ** 2
            + 10.1 * ((y - 1) ** 2 + (w - 1) ** 2) + 19.8 * (y - 1) * (w - 1))

# ---- 07 Trid n, f* = -50 at n=6 ----
def _trid(v):
    return float(np.sum((v - 1) ** 2) - np.sum(v[1:] * v[:-1]))

# ---- 08 bound-active quadratic on [0,1]^2, f* = 5 ----
def _bq(v):   x, y = v; return (x - 2) ** 2 + (y - 3) ** 2

# ---- 09 non-smooth L1 sum, f* = 0 at (1..n) ----
def _l1(v):   return float(np.sum(np.abs(v - (1.0 + np.arange(v.size)))))

# ---- separable / ill-conditioned diagonal quadratics ----
def _sphere(v):  return float(np.sum(v * v))
def _ill(n, c):
    w = 10.0 ** (c * np.arange(n) / (n - 1))
    return (lambda v: float(np.sum(w * v * v)), np.ones(n))


ill10c4f, ill10c4x = _ill(10, 4)
ill10c6f, ill10c6x = _ill(10, 6)

bench("01 Rosenbrock2D", lambda: _solve(_rb, [-1.2, 1.0]))
check("01 Rosenbrock2D", int(np.floor(1e6 * _solve(_rb, [-1.2, 1.0]).fun + 0.5)))

bench("02 Beale", lambda: _solve(_beale, [1.0, 1.0]))
check("02 Beale", int(np.floor(1e6 * _solve(_beale, [1.0, 1.0]).fun + 0.5)))

bench("03 Booth", lambda: _solve(_booth, [0.0, 0.0]))
check("03 Booth", int(np.floor(1e6 * _solve(_booth, [0.0, 0.0]).fun + 0.5)))

bench("04 Matyas", lambda: _solve(_mat, [3.0, 3.0]))
check("04 Matyas", int(np.floor(1e6 * _solve(_mat, [3.0, 3.0]).fun + 0.5)))

bench("05 Powell-singular", lambda: _solve(_pow, [3.0, -1.0, 0.0, 1.0]))
check("05 Powell-singular", int(np.floor(1e6 * _solve(_pow, [3.0, -1.0, 0.0, 1.0]).fun + 0.5)))

bench("06 Wood", lambda: _solve(_wood, [-3.0, -1.0, -3.0, -1.0]))
check("06 Wood", int(np.floor(1e6 * _solve(_wood, [-3.0, -1.0, -3.0, -1.0]).fun + 0.5)))

bench("07 Trid n6", lambda: _solve(_trid, np.zeros(6)))
check("07 Trid n6", int(np.floor(_solve(_trid, np.zeros(6)).fun + 0.5)))

bench("08 Bound-corner", lambda: _solve(_bq, [0.0, 0.0], bounds=[(0.0, 1.0), (0.0, 1.0)]))
check("08 Bound-corner", int(np.floor(_solve(_bq, [0.0, 0.0], bounds=[(0.0, 1.0), (0.0, 1.0)]).fun + 0.5)))

bench("09 Nonsmooth-L1 n5", lambda: _solve(_l1, np.zeros(5)))
check("09 Nonsmooth-L1 n5", int(np.floor(1e6 * _solve(_l1, np.zeros(5)).fun + 0.5)))

bench("10 Sphere n10", lambda: _solve(_sphere, np.ones(10)))
check("10 Sphere n10", int(np.floor(1e6 * _solve(_sphere, np.ones(10)).fun + 0.5)))

bench("11 Sphere n20", lambda: _solve(_sphere, np.ones(20)))
check("11 Sphere n20", int(np.floor(1e6 * _solve(_sphere, np.ones(20)).fun + 0.5)))

bench("12 illcond n10 c1e4", lambda: _solve(ill10c4f, ill10c4x))
check("12 illcond n10 c1e4", int(np.floor(1e6 * _solve(ill10c4f, ill10c4x).fun + 0.5)))

bench("13 illcond n10 c1e6", lambda: _solve(ill10c6f, ill10c6x))
check("13 illcond n10 c1e6", int(np.floor(1e6 * _solve(ill10c6f, ill10c6x).fun + 0.5)))
