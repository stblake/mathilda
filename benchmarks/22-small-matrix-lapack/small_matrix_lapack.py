#!/usr/bin/env python3
"""Experiment 22 -- Small matrices: does LAPACK get reached at all? (numpy column)

Same six kernels as ``small_matrix_lapack.m``, same order and repetition counts.

ROADMAP ITEM 3.  A small matrix sits below the packing threshold, so it is a
plain nested List, so a dispatch that looks for an *existing buffer* misses and
the generic symbolic path runs.  These rows are small-matrix operations repeated
many times -- the shape a Kalman filter or a ray tracer actually produces, and
the shape a per-value threshold cannot see.

numpy has the same per-call overhead problem in miniature (a 6x6 ``inv`` is
dominated by Python and LAPACK call setup, not by the 216 flops), which makes
this a fair overhead-vs-overhead comparison rather than a flops comparison.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import numpy as np

from harness import bench, check, require
from data import spd_int_matrix

require(["numpy.linalg:inv", "numpy.linalg:det", "numpy.linalg:solve",
         "numpy.linalg:eigvals"])

m6 = spd_int_matrix(6)
m3 = spd_int_matrix(3)
v6 = np.arange(1.0, 7.0)


def r(x, p):
    return int(math.floor(float(x) * 10 ** p + 0.5))


bench("Inverse 6x6 x 2000",
      lambda: [np.linalg.inv(m6) for _ in range(2000)])
check("Inverse 6x6 x 2000", r(np.linalg.inv(m6).sum(), 6))

bench("Det 6x6 x 5000", lambda: [np.linalg.det(m6) for _ in range(5000)])
check("Det 6x6 x 5000", r(np.linalg.det(m6), 4))

bench("LinearSolve 6x6 x 2000",
      lambda: [np.linalg.solve(m6, v6) for _ in range(2000)])
check("LinearSolve 6x6 x 2000", r(np.linalg.solve(m6, v6).sum(), 6))

bench("Dot 6x6 x 6x6 x 10000", lambda: [m6 @ m6 for _ in range(10000)])
check("Dot 6x6 x 6x6 x 10000", r((m6 @ m6).sum(), 4))

bench("Inverse 3x3 x 5000",
      lambda: [np.linalg.inv(m3) for _ in range(5000)])
check("Inverse 3x3 x 5000", r(np.linalg.inv(m3).sum(), 6))

bench("Eigenvalues 6x6 x 500",
      lambda: [np.linalg.eigvals(m6) for _ in range(500)])
check("Eigenvalues 6x6 x 500", r(np.sort(np.linalg.eigvals(m6).real).sum(), 4))
