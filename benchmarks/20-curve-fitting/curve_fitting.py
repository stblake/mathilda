#!/usr/bin/env python3
"""Experiment 20 -- Curve fitting and regression (numpy/scipy column).

Same five kernels as ``curve_fitting.m``, same order and sizes.

Separate from experiment 15 (optimisation) because linear least squares is a
DIRECT solve, not an iteration -- different code, different failure mode.  All
data is deterministic, so the three systems fit the same points.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import numpy as np

from harness import bench, check, require, seed
from data import spd_int_matrix

require(["numpy:polyfit", "numpy.linalg:lstsq", "scipy.linalg:lstsq"])


def r(x, p):
    return int(math.floor(float(x) * 10 ** p + 0.5))


t1 = np.arange(1, 5001) / 100.0
y1 = 3.0 * t1 + 2.0
bench("Fit linear, 5000 points", lambda: np.polyfit(t1, y1, 1))
check("Fit linear, 5000 points", r(np.polyval(np.polyfit(t1, y1, 1), 7.0), 4))

t2 = np.arange(1, 5001) / 500.0
y2 = t2 ** 2 - t2 + 1.0
bench("Fit quadratic, 5000 points", lambda: np.polyfit(t2, y2, 2))
check("Fit quadratic, 5000 points",
      r(np.polyval(np.polyfit(t2, y2, 2), 3.0), 4))

bench("Fit degree-8 polynomial, 5000 points",
      lambda: np.polyfit(t2, y2, 8), reps=1)
check("Fit degree-8 polynomial, 5000 points",
      r(np.polyval(np.polyfit(t2, y2, 8), 3.0), 2))

seed()
m = np.random.random((2000, 20))
rhs = np.random.random(2000)
bench("LeastSquares 2000x20", lambda: np.linalg.lstsq(m, rhs, rcond=None))
check("LeastSquares 2000x20",
      r(np.linalg.lstsq(spd_int_matrix(8), np.arange(1.0, 9.0),
                        rcond=None)[0].sum(), 4))


def _trigfit():
    A = np.column_stack([np.ones_like(t1), np.sin(t1), np.cos(t1)])
    return np.linalg.lstsq(A, y1, rcond=None)[0]


bench("Fit trig basis, 5000 points", _trigfit)
_c = _trigfit()
check("Fit trig basis, 5000 points",
      r(_c[0] + _c[1] * math.sin(1.0) + _c[2] * math.cos(1.0), 2))
