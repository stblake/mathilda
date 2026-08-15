#!/usr/bin/env python3
"""Experiment 69 -- COBYLA (derivative-free constrained) optimization (scipy column).

Same 4 problems as ``cobyla_constrained.m``, same labels, same order, solved with
scipy.optimize.minimize(method="COBYLA"). Both solvers are derivative-free (no
Jacobian either side), so this is a fair derivative-free-vs-derivative-free race.
scipy's COBYLA is inequality-only (fun(x) >= 0 feasible); bounds are passed via
bounds= (scipy >= 1.11). The check is Round[10^3 * f] at the optimum; both must
agree or the join CHECK-FAILs.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import minimize

from harness import bench, check, require

require(["scipy.optimize:minimize"])

_OPTS = {"maxiter": 3000, "tol": 1e-8}


# ---- 01 QP tutorial: min (x-1)^2+(y-2.5)^2, 3 linear inequalities + bounds
def _qp(v):   x, y = v; return (x - 1) ** 2 + (y - 2.5) ** 2
_qp_cons = [
    {"type": "ineq", "fun": lambda v: v[0] - 2 * v[1] + 2},
    {"type": "ineq", "fun": lambda v: -v[0] - 2 * v[1] + 6},
    {"type": "ineq", "fun": lambda v: -v[0] + 2 * v[1] + 2},
]

def _solve_qp():
    return minimize(_qp, [2.0, 0.0], method="COBYLA",
                    bounds=[(0.0, None), (0.0, None)], constraints=_qp_cons, options=_OPTS)

bench("01 QP tutorial", _solve_qp)
check("01 QP tutorial", int(np.floor(1e3 * _solve_qp().fun + 0.5)))


# ---- 02 Corner LP: min x+y s.t. 3x+2y>=7, x>=0, y>=0 -> 7/3 at (7/3,0)
def _corner(v):  return v[0] + v[1]
_corner_cons = [{"type": "ineq", "fun": lambda v: 3 * v[0] + 2 * v[1] - 7}]

def _solve_corner():
    return minimize(_corner, [1.0, 1.0], method="COBYLA",
                    bounds=[(0.0, None), (0.0, None)], constraints=_corner_cons, options=_OPTS)

bench("02 Corner LP", _solve_corner)
check("02 Corner LP", int(np.floor(1e3 * _solve_corner().fun + 0.5)))


# ---- 03 Nonlinear inequality: min x^2+y^2 s.t. x y >= 1 -> 2 at (1,1)
def _nl(v):   return v[0] ** 2 + v[1] ** 2
_nl_cons = [{"type": "ineq", "fun": lambda v: v[0] * v[1] - 1.0}]

def _solve_nl():
    return minimize(_nl, [2.0, 0.5], method="COBYLA", constraints=_nl_cons, options=_OPTS)

bench("03 Nonlinear-ineq xy>=1", _solve_nl)
check("03 Nonlinear-ineq xy>=1", int(np.floor(1e3 * _solve_nl().fun + 0.5)))


# ---- 04 Constrained quadratic n=5: min sum (x_i - i)^2 s.t. sum x <= 10 -> 5
def _make_cq(n):
    idx = np.arange(1, n + 1, dtype=float)
    fun = lambda v: float(np.sum((np.asarray(v) - idx) ** 2))
    cons = [{"type": "ineq", "fun": lambda v: 10.0 - float(np.sum(v))}]
    x0 = np.zeros(n)

    def solve():
        return minimize(fun, x0, method="COBYLA", constraints=cons, options=_OPTS)
    return solve

_cq5 = _make_cq(5)
bench("04 Constrained-quad n5", _cq5)
check("04 Constrained-quad n5", int(np.floor(1e3 * _cq5().fun + 0.5)))
