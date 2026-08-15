#!/usr/bin/env python3
"""Experiment 68 -- SLSQP constrained optimization (scipy column).

Same 6 problems as ``slsqp_constrained.m``, same labels, same order, solved with
scipy.optimize.minimize(method="SLSQP") and analytic Jacobians for both the
objective and every constraint (so the race is solver-vs-solver -- both are
gradient-based SQP). scipy's convention: inequality constraints are fun(x) >= 0.
The check is the objective at the optimum rounded to 6 digits (10^3 for HS71,
whose optimum is irrational): both systems must reach it or the join CHECK-FAILs.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import minimize

from harness import bench, check, require

require(["scipy.optimize:minimize"])

_OPTS = {"maxiter": 500, "ftol": 1e-12}


# ---- 01 QP tutorial: min (x-1)^2+(y-2.5)^2, five linear inequalities + bounds
def _qp(v):   x, y = v; return (x - 1) ** 2 + (y - 2.5) ** 2
def _qpg(v):  x, y = v; return np.array([2 * (x - 1), 2 * (y - 2.5)])
_qp_cons = [
    {"type": "ineq", "fun": lambda v: v[0] - 2 * v[1] + 2,  "jac": lambda v: np.array([1.0, -2.0])},
    {"type": "ineq", "fun": lambda v: -v[0] - 2 * v[1] + 6, "jac": lambda v: np.array([-1.0, -2.0])},
    {"type": "ineq", "fun": lambda v: -v[0] + 2 * v[1] + 2, "jac": lambda v: np.array([-1.0, 2.0])},
]

def _solve_qp():
    return minimize(_qp, [2.0, 0.0], jac=_qpg, method="SLSQP",
                    bounds=[(0.0, None), (0.0, None)], constraints=_qp_cons, options=_OPTS)

bench("01 QP tutorial", _solve_qp)
check("01 QP tutorial", int(np.floor(1e6 * _solve_qp().fun + 0.5)))


# ---- 02 HS71: min x1 x4 (x1+x2+x3) + x3; ineq x1x2x3x4>=25; eq sum sq ==40; 1<=xi<=5
def _hs71(v):
    x1, x2, x3, x4 = v
    return x1 * x4 * (x1 + x2 + x3) + x3
def _hs71g(v):
    x1, x2, x3, x4 = v
    return np.array([x4 * (2 * x1 + x2 + x3),
                     x1 * x4,
                     x1 * x4 + 1.0,
                     x1 * (x1 + x2 + x3)])
_hs71_cons = [
    {"type": "ineq",
     "fun": lambda v: v[0] * v[1] * v[2] * v[3] - 25.0,
     "jac": lambda v: np.array([v[1] * v[2] * v[3], v[0] * v[2] * v[3],
                                v[0] * v[1] * v[3], v[0] * v[1] * v[2]])},
    {"type": "eq",
     "fun": lambda v: v[0] ** 2 + v[1] ** 2 + v[2] ** 2 + v[3] ** 2 - 40.0,
     "jac": lambda v: 2.0 * np.asarray(v)},
]

def _solve_hs71():
    return minimize(_hs71, [1.0, 5.0, 5.0, 1.0], jac=_hs71g, method="SLSQP",
                    bounds=[(1.0, 5.0)] * 4, constraints=_hs71_cons, options=_OPTS)

bench("02 HS71", _solve_hs71)
check("02 HS71", int(np.floor(1e3 * _solve_hs71().fun + 0.5)))


# ---- 03 nonlinear equality: min x^2+y^2 s.t. x y == 1 -> f* = 2 at (1,1)
def _nl(v):   x, y = v; return x * x + y * y
def _nlg(v):  x, y = v; return np.array([2 * x, 2 * y])
_nl_cons = [{"type": "eq", "fun": lambda v: v[0] * v[1] - 1.0,
             "jac": lambda v: np.array([v[1], v[0]])}]

def _solve_nl():
    return minimize(_nl, [2.0, 0.5], jac=_nlg, method="SLSQP",
                    constraints=_nl_cons, options=_OPTS)

bench("03 Nonlinear-eq xy=1", _solve_nl)
check("03 Nonlinear-eq xy=1", int(np.floor(1e6 * _solve_nl().fun + 0.5)))


# ---- 04-06 equal-weight simplex QP: min sum w^2 s.t. sum w == 1, w >= 0 -> 1/n
def _make_simplex(n):
    fun = lambda v: float(np.dot(v, v))
    jac = lambda v: 2.0 * np.asarray(v)
    cons = [{"type": "eq", "fun": lambda v: float(np.sum(v)) - 1.0,
             "jac": lambda v: np.ones(n)}]
    x0 = np.full(n, 1.0 / n)
    bounds = [(0.0, None)] * n

    def solve():
        return minimize(fun, x0, jac=jac, method="SLSQP",
                        bounds=bounds, constraints=cons, options=_OPTS)
    return solve

_s10 = _make_simplex(10)
_s50 = _make_simplex(50)
_s100 = _make_simplex(100)

bench("04 Simplex QP n10", _s10)
check("04 Simplex QP n10", int(np.floor(1e6 * _s10().fun + 0.5)))

bench("05 Simplex QP n50", _s50)
check("05 Simplex QP n50", int(np.floor(1e6 * _s50().fun + 0.5)))

bench("06 Simplex QP n100", _s100)
check("06 Simplex QP n100", int(np.floor(1e6 * _s100().fun + 0.5)))
