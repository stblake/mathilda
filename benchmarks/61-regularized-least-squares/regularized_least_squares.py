#!/usr/bin/env python3
"""Experiment 61 -- Least squares & regularized fitting (numpy/scipy column).

Same kernels as ``regularized_least_squares.m``, same order and labels.

Baseline is compiled: numpy.linalg.lstsq / scipy.linalg.lstsq (LAPACK gelsd),
numpy.polyfit, and a numpy closed-form ridge (X^T X + lambda I)^-1 X^T y, which
Mathilda's Tikhonov FitRegularization matches.  FindFit is absent in Mathilda,
so scipy.optimize.curve_fit stands in for it (the .m side SKIPs it -> ABSENT).

Timing on large random systems; checks on the small deterministic Vandermonde.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math
import numpy as np
import scipy.linalg as sla
from scipy.optimize import curve_fit

from harness import bench, check, require, seed

require(["numpy.linalg:lstsq", "scipy.linalg:lstsq", "numpy:polyfit",
         "scipy.optimize:curve_fit"])

seed()
a1 = np.random.random((2000, 50)); b1 = np.random.random(2000)
a2 = np.random.random((8000, 200)); b2 = np.random.random(8000)
xbig = np.arange(1, 2001) / 100.0
ybig = 2 + 3 * xbig + xbig ** 2 + np.sin(np.arange(1, 2001, dtype=float))

# deterministic check system
xc = np.arange(1, 13, dtype=float)
des = np.array([[x ** j for j in range(5)] for x in xc])
bc = np.array([2 + 3 * x + x * x + (int(x) % 3) for x in xc], dtype=float)


def r(x, p):
    return int(math.floor(float(x) * 10 ** p + 0.5))


def ls_check():
    return np.linalg.lstsq(des, bc, rcond=None)[0].sum()


bench("LeastSquares overdetermined 2000x50", lambda: np.linalg.lstsq(a1, b1, rcond=None))
check("LeastSquares overdetermined 2000x50", r(ls_check(), 6))

bench("LeastSquares overdetermined 8000x200", lambda: sla.lstsq(a2, b2), reps=1)
check("LeastSquares overdetermined 8000x200", r(ls_check(), 6))

bench("Fit polynomial degree 4", lambda: np.polyfit(xbig, ybig, 4))
check("Fit polynomial degree 4", r(np.polyval(np.polyfit(xc, bc, 4), 7.0), 6))


def ridge(X, y, lam):
    return np.linalg.solve(X.T @ X + lam * np.eye(X.shape[1]), X.T @ y)


Xbig = np.array([xbig ** j for j in range(5)]).T
bench("Fit ridge (Tikhonov) degree 4", lambda: ridge(Xbig, ybig, 0.5))
_a = ridge(des, bc, 0.5)
check("Fit ridge (Tikhonov) degree 4", r(sum(_a[j] * 7.0 ** j for j in range(5)), 6))

# FindFit is ABSENT in Mathilda; curve_fit measures the Python side.
dataC_x = xc
dataC_y = bc
bench("FindFit exponential model",
      lambda: curve_fit(lambda x, a, b: a * np.exp(b * x), dataC_x, dataC_y, p0=[1.0, 0.1]))
check("FindFit exponential model",
      r(curve_fit(lambda x, a, b: a * np.exp(b * x), dataC_x, dataC_y, p0=[1.0, 0.1])[0][0], 3))
