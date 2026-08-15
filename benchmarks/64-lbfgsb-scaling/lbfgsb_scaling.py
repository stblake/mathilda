#!/usr/bin/env python3
"""Experiment 64 -- L-BFGS-B scaling (scipy.optimize.minimize column).

Same cases as ``lbfgsb_scaling.m``, same labels, same order. The check is the
objective value at the optimum, rounded to an integer -- invariant under the
descent path, so the two systems join cleanly. scipy is given the ANALYTIC
Jacobian (jac=) to match Mathilda's exact compiled gradient: solver-vs-solver.

The matching scipy method mirrors the Mathilda method: "LBFGSB" races
method="L-BFGS-B", "QuasiNewton" races the full-memory method="BFGS".
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import minimize

from harness import bench, check, require

require(["scipy.optimize:minimize"])

_LB_OPTS = {"maxiter": 15000, "maxfun": 200000, "ftol": 1e-15, "gtol": 1e-8}
_BFGS_OPTS = {"maxiter": 15000, "gtol": 1e-8}


# ---- ill-conditioned diagonal quadratic (condition 1e2), min 0 at origin ----
def _cq_w(n):
    return 10.0 ** (2.0 * np.arange(n) / (n - 1))


def _solve_cq(n, method):
    w = _cq_w(n)
    fun = lambda x: float(np.sum(w * x * x))
    jac = lambda x: 2.0 * w * x
    x0 = np.ones(n)
    opts = _LB_OPTS if method == "L-BFGS-B" else _BFGS_OPTS
    return minimize(fun, x0, jac=jac, method=method, options=opts)


# ---- separable box quadratic on [0,1]^n, optimum the corner (all 1), value n --
def _solve_box(n):
    fun = lambda x: float(np.sum((x - 2.0) ** 2))
    jac = lambda x: 2.0 * (x - 2.0)
    x0 = np.zeros(n)
    return minimize(fun, x0, jac=jac, method="L-BFGS-B",
                    bounds=[(0.0, 1.0)] * n, options=_LB_OPTS)


# ---- Rosenbrock (extended sum form), global optimum 0 at all-ones ----------
def _rosen(x):
    return float(np.sum(100.0 * (x[1:] - x[:-1] ** 2) ** 2 + (1.0 - x[:-1]) ** 2))


def _rosen_grad(x):
    n = x.size
    g = np.zeros(n)
    xm = x[:-1]
    g[:-1] = -400.0 * xm * (x[1:] - xm ** 2) - 2.0 * (1.0 - xm)
    g[1:] += 200.0 * (x[1:] - xm ** 2)
    return g


def _solve_rosen(x0):
    return minimize(_rosen, x0, jac=_rosen_grad, method="L-BFGS-B", options=_LB_OPTS)


# ---- C-series: L-BFGS-B, scaling in n --------------------------------------
bench("C50 quadratic LBFGSB n=50", lambda: _solve_cq(50, "L-BFGS-B"))
check("C50 quadratic LBFGSB n=50", int(np.floor(_solve_cq(50, "L-BFGS-B").fun + 0.5)))

bench("C200 quadratic LBFGSB n=200", lambda: _solve_cq(200, "L-BFGS-B"))
check("C200 quadratic LBFGSB n=200", int(np.floor(_solve_cq(200, "L-BFGS-B").fun + 0.5)))

bench("C1000 quadratic LBFGSB n=1000", lambda: _solve_cq(1000, "L-BFGS-B"))
check("C1000 quadratic LBFGSB n=1000", int(np.floor(_solve_cq(1000, "L-BFGS-B").fun + 0.5)))

# ---- Q-series: full-memory BFGS on the same quadratic ----------------------
bench("Q50 quadratic QuasiNewton n=50", lambda: _solve_cq(50, "BFGS"))
check("Q50 quadratic QuasiNewton n=50", int(np.floor(_solve_cq(50, "BFGS").fun + 0.5)))

bench("Q200 quadratic QuasiNewton n=200", lambda: _solve_cq(200, "BFGS"))
check("Q200 quadratic QuasiNewton n=200", int(np.floor(_solve_cq(200, "BFGS").fun + 0.5)))

bench("Q1000 quadratic QuasiNewton n=1000", lambda: _solve_cq(1000, "BFGS"))
check("Q1000 quadratic QuasiNewton n=1000", int(np.floor(_solve_cq(1000, "BFGS").fun + 0.5)))

# ---- B-series: bound-constrained box (the "-B" path) -----------------------
bench("B200 bound-box LBFGSB n=200", lambda: _solve_box(200))
check("B200 bound-box LBFGSB n=200", int(np.floor(_solve_box(200).fun + 0.5)))

bench("B1000 bound-box LBFGSB n=1000", lambda: _solve_box(1000))
check("B1000 bound-box LBFGSB n=1000", int(np.floor(_solve_box(1000).fun + 0.5)))

# ---- R-series: small-n Rosenbrock ------------------------------------------
bench("R2 rosenbrock LBFGSB n=2", lambda: _solve_rosen(np.array([-1.2, 1.0])))
check("R2 rosenbrock LBFGSB n=2", int(np.floor(_solve_rosen(np.array([-1.2, 1.0])).fun + 0.5)))

bench("R5 rosenbrock LBFGSB n=5", lambda: _solve_rosen(np.full(5, 0.5)))
check("R5 rosenbrock LBFGSB n=5", int(np.floor(_solve_rosen(np.full(5, 0.5)).fun + 0.5)))
