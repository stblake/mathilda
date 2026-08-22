#!/usr/bin/env python3
"""Experiment 72 -- dogleg unconstrained optimization (scipy column).

Same 6 problems as ``dogleg.m``, same labels, same order, solved with
scipy.optimize.minimize(method="dogleg") and the analytic Jacobian AND Hessian,
so the race is solver-vs-solver (both Powell dogleg trust region).

dogleg needs a POSITIVE-DEFINITE model Hessian (scipy raises LinAlgError
otherwise), so the problem set is PD throughout: ill-conditioned diagonal
quadratics (condition 1, 1e2, 1e4, 1e6) and the Trid function (constant SPD
Hessian). The check is the objective at the optimum rounded to 10^6 (0 for the
f*=0 quadratics, the exact integer for Trid).
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import minimize

from harness import bench, check, require

require(["scipy.optimize:minimize"])

_METHOD = "dogleg"
_OPTS = {"maxiter": 5000, "gtol": 1e-10}


def _solve(fun, jac, hess, x0):
    return minimize(fun, np.asarray(x0, float), jac=jac, hess=hess,
                    method=_METHOD, options=_OPTS)


# ---- ill-conditioned diagonal quadratic  f = sum 10^(c(i-1)/(n-1)) x_i^2 ----
def _quad_case(n, c):
    a = np.array([10.0 ** (c * i / (n - 1)) for i in range(n)])
    fun = lambda v: float(np.dot(a, v * v))
    jac = lambda v: 2.0 * a * v
    hess = lambda v: np.diag(2.0 * a)
    x0 = np.ones(n)
    return lambda: _solve(fun, jac, hess, x0)


# ---- Trid: f = sum (x_i-1)^2 - sum x_i x_{i-1}; min -n(n+4)(n-1)/6 ----
def _trid_case(n):
    def fun(v):
        return float(np.sum((v - 1.0) ** 2) - np.sum(v[1:] * v[:-1]))
    def jac(v):
        g = 2.0 * (v - 1.0)
        g[:-1] -= v[1:]
        g[1:] -= v[:-1]
        return g
    def hess(v):
        H = 2.0 * np.eye(n)
        for i in range(n - 1):
            H[i, i + 1] = -1.0
            H[i + 1, i] = -1.0
        return H
    x0 = np.zeros(n)
    return lambda: _solve(fun, jac, hess, x0)


_q0, _q2, _q4, _q6 = _quad_case(10, 0), _quad_case(10, 2), _quad_case(10, 4), _quad_case(10, 6)
_t10, _t20 = _trid_case(10), _trid_case(20)

bench("01 Sphere n10", _q0)
check("01 Sphere n10", int(np.floor(1e6 * _q0().fun + 0.5)))

bench("02 Illcond quad n10 c1e2", _q2)
check("02 Illcond quad n10 c1e2", int(np.floor(1e6 * _q2().fun + 0.5)))

bench("03 Illcond quad n10 c1e4", _q4)
check("03 Illcond quad n10 c1e4", int(np.floor(1e6 * _q4().fun + 0.5)))

bench("04 Illcond quad n10 c1e6", _q6)
check("04 Illcond quad n10 c1e6", int(np.floor(1e6 * _q6().fun + 0.5)))

bench("05 Trid n10", _t10)
check("05 Trid n10", int(np.floor(1e6 * _t10().fun + 0.5)))

bench("06 Trid n20", _t20)
check("06 Trid n20", int(np.floor(1e6 * _t20().fun + 0.5)))
