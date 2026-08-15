#!/usr/bin/env python3
"""Experiment 70 -- COBYQA (derivative-free quadratic, constrained) (scipy column).

Same 4 problems as ``cobyqa_constrained.m``, same labels, same order, solved with
scipy.optimize.minimize(method="COBYQA") (Ragonneau & Zhang 2023). Both solvers
are derivative-free (no Jacobian either side). Unlike COBYLA, COBYQA handles
equality constraints natively, so this suite mixes eq + ineq + bounds. The check
is Round[10^3 * f] at the optimum (10^2 for HS71, whose optimum is irrational);
both must agree or the join CHECK-FAILs.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import minimize

from harness import bench, check, require

require(["scipy.optimize:minimize"])

_OPTS = {"maxiter": 3000}


# ---- 01 HS71: min x1 x4 (x1+x2+x3) + x3; ineq x1x2x3x4>=25; eq sumsq==40; 1<=xi<=5
def _hs71(v):
    x1, x2, x3, x4 = v
    return x1 * x4 * (x1 + x2 + x3) + x3
_hs71_cons = [
    {"type": "ineq", "fun": lambda v: v[0] * v[1] * v[2] * v[3] - 25.0},
    {"type": "eq",   "fun": lambda v: v[0] ** 2 + v[1] ** 2 + v[2] ** 2 + v[3] ** 2 - 40.0},
]

def _solve_hs71():
    return minimize(_hs71, [1.0, 5.0, 5.0, 1.0], method="COBYQA",
                    bounds=[(1.0, 5.0)] * 4, constraints=_hs71_cons, options=_OPTS)

bench("01 HS71", _solve_hs71)
check("01 HS71", int(np.floor(1e2 * _solve_hs71().fun + 0.5)))


# ---- 02 nonlinear equality: min x^2+y^2 s.t. x y == 1 -> 2 at (1,1)
def _nl(v):   return v[0] ** 2 + v[1] ** 2
_nl_cons = [{"type": "eq", "fun": lambda v: v[0] * v[1] - 1.0}]

def _solve_nl():
    return minimize(_nl, [1.5, 0.8], method="COBYQA", constraints=_nl_cons, options=_OPTS)

bench("02 Nonlinear-eq xy=1", _solve_nl)
check("02 Nonlinear-eq xy=1", int(np.floor(1e3 * _solve_nl().fun + 0.5)))


# ---- 03 constrained Rosenbrock: min (1-x)^2+100(y-x^2)^2 s.t. x+y<=1
def _cros(v):  x, y = v; return (1 - x) ** 2 + 100 * (y - x * x) ** 2
_cros_cons = [{"type": "ineq", "fun": lambda v: 1.0 - v[0] - v[1]}]

def _solve_cros():
    return minimize(_cros, [-1.0, 1.0], method="COBYQA", constraints=_cros_cons, options=_OPTS)

bench("03 Constrained-Rosenbrock", _solve_cros)
check("03 Constrained-Rosenbrock", int(np.floor(1e3 * _solve_cros().fun + 0.5)))


# ---- 04 equal-weight simplex QP n=8: min sum w^2 s.t. sum w == 1, w >= 0 -> 1/8
def _make_simplex(n):
    fun = lambda v: float(np.dot(v, v))
    cons = [{"type": "eq", "fun": lambda v: float(np.sum(v)) - 1.0}]
    x0 = np.full(n, 1.0 / n)
    bounds = [(0.0, None)] * n

    def solve():
        return minimize(fun, x0, method="COBYQA", bounds=bounds, constraints=cons, options=_OPTS)
    return solve

_s8 = _make_simplex(8)
bench("04 Simplex QP n8", _s8)
check("04 Simplex QP n8", int(np.floor(1e3 * _s8().fun + 0.5)))
