#!/usr/bin/env python3
"""Experiment 67 -- TNC (truncated Newton) scaling + accuracy (scipy column).

Same 11 test functions as ``tnc_scaling.m``, same labels, same order, solved with
scipy.optimize.minimize(method="TNC") and the analytic Jacobian (jac=), so the
race is solver-vs-solver (both gradient-based truncated Newton). The check is the
objective at the optimum rounded to 6 digits (0 for the f*=0 problems) or the
exact integer optimum (Trid -50, bound corner 5, many-active 160).
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import minimize

from harness import bench, check, require

require(["scipy.optimize:minimize"])

_OPTS = {"maxfun": 200000, "ftol": 1e-14, "xtol": 1e-12, "gtol": 1e-10}


def _solve(fun, jac, x0, bounds=None):
    return minimize(fun, np.asarray(x0, float), jac=jac, method="TNC",
                    bounds=bounds, options=_OPTS)


# ---- 01 Rosenbrock 2D ----
def _rb(v):   x, y = v; return (1 - x) ** 2 + 100 * (y - x * x) ** 2
def _rbg(v):  x, y = v; return np.array([-2 * (1 - x) - 400 * x * (y - x * x), 200 * (y - x * x)])

# ---- extended Rosenbrock n (contiguous coupling, order z1..zn) ----
def _erb(v):
    return float(np.sum(100 * (v[1:] - v[:-1] ** 2) ** 2 + (1 - v[:-1]) ** 2))
def _erbg(v):
    g = np.zeros_like(v); x0 = v[:-1]; x1 = v[1:]
    g[:-1] += -400 * x0 * (x1 - x0 ** 2) - 2 * (1 - x0)
    g[1:]  += 200 * (x1 - x0 ** 2)
    return g

# ---- ill-conditioned diagonal quadratic ----
def _ill(n, c):
    w = 10.0 ** (c * np.arange(n) / (n - 1))
    return (lambda v: float(np.sum(w * v * v)), lambda v: 2.0 * w * v, np.ones(n))

# ---- Trid n ----
def _trid(v):
    return float(np.sum((v - 1) ** 2) - np.sum(v[1:] * v[:-1]))
def _tridg(v):
    g = 2 * (v - 1); g[1:] -= v[:-1]; g[:-1] -= v[1:]; return g

# ---- bound-active quadratic on [0,1]^2 ----
def _bq(v):   x, y = v; return (x - 2) ** 2 + (y - 3) ** 2
def _bqg(v):  x, y = v; return np.array([2 * (x - 2), 2 * (y - 3)])

# ---- many-active: sum (z_i-5)^2, z_i <= 1 ----
def _ma(v):   return float(np.sum((v - 5) ** 2))
def _mag(v):  return 2.0 * (v - 5)

# ---- Booth ----
def _booth(v):  x, y = v; return (x + 2 * y - 7) ** 2 + (2 * x + y - 5) ** 2
def _boothg(v):
    x, y = v
    return np.array([2 * (x + 2 * y - 7) + 4 * (2 * x + y - 5),
                     4 * (x + 2 * y - 7) + 2 * (2 * x + y - 5)])

# ---- Wood 4D (order x, y, z, w) ----
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


ill10c6f, ill10c6g, ill10c6x = _ill(10, 6)
ill50c4f, ill50c4g, ill50c4x = _ill(50, 4)
ill100c4f, ill100c4g, ill100c4x = _ill(100, 4)

bench("01 Rosenbrock2D", lambda: _solve(_rb, _rbg, [-1.2, 1.0]))
check("01 Rosenbrock2D", int(np.floor(1e6 * _solve(_rb, _rbg, [-1.2, 1.0]).fun + 0.5)))

bench("02 Rosenbrock n10", lambda: _solve(_erb, _erbg, -1.2 * np.ones(10)))
check("02 Rosenbrock n10", int(np.floor(1e6 * _solve(_erb, _erbg, -1.2 * np.ones(10)).fun + 0.5)))

bench("03 Rosenbrock n30", lambda: _solve(_erb, _erbg, -1.2 * np.ones(30)))
check("03 Rosenbrock n30", int(np.floor(1e6 * _solve(_erb, _erbg, -1.2 * np.ones(30)).fun + 0.5)))

bench("04 illcond n10 c1e6", lambda: _solve(ill10c6f, ill10c6g, ill10c6x))
check("04 illcond n10 c1e6", int(np.floor(1e6 * _solve(ill10c6f, ill10c6g, ill10c6x).fun + 0.5)))

bench("05 illcond n50 c1e4", lambda: _solve(ill50c4f, ill50c4g, ill50c4x))
check("05 illcond n50 c1e4", int(np.floor(1e6 * _solve(ill50c4f, ill50c4g, ill50c4x).fun + 0.5)))

bench("06 illcond n100 c1e4", lambda: _solve(ill100c4f, ill100c4g, ill100c4x))
check("06 illcond n100 c1e4", int(np.floor(1e6 * _solve(ill100c4f, ill100c4g, ill100c4x).fun + 0.5)))

bench("07 Trid n6", lambda: _solve(_trid, _tridg, np.zeros(6)))
check("07 Trid n6", int(np.floor(_solve(_trid, _tridg, np.zeros(6)).fun + 0.5)))

bench("08 Bound-corner", lambda: _solve(_bq, _bqg, [0.0, 0.0], bounds=[(0.0, 1.0), (0.0, 1.0)]))
check("08 Bound-corner", int(np.floor(_solve(_bq, _bqg, [0.0, 0.0], bounds=[(0.0, 1.0), (0.0, 1.0)]).fun + 0.5)))

bench("09 Many-active n10", lambda: _solve(_ma, _mag, np.zeros(10), bounds=[(-50.0, 1.0)] * 10))
check("09 Many-active n10", int(np.floor(_solve(_ma, _mag, np.zeros(10), bounds=[(-50.0, 1.0)] * 10).fun + 0.5)))

bench("10 Booth", lambda: _solve(_booth, _boothg, [0.0, 0.0]))
check("10 Booth", int(np.floor(1e6 * _solve(_booth, _boothg, [0.0, 0.0]).fun + 0.5)))

bench("11 Wood", lambda: _solve(_wood, _woodg, [-3.0, -1.0, -3.0, -1.0]))
check("11 Wood", int(np.floor(1e6 * _solve(_wood, _woodg, [-3.0, -1.0, -3.0, -1.0]).fun + 0.5)))
