#!/usr/bin/env python3
"""Experiment 18 -- Dense linear algebra drivers (scipy.linalg column).

Same seven kernels as ``dense_linalg_drivers.m``, same order and sizes.

On this host numpy links the SAME Accelerate BLAS Mathilda does, so on these rows
all three systems call byte-identical kernels and any spread is PURE OVERHEAD --
the cleanest measurement in the suite.

Checks use the exactly reproducible SPD integer matrix from data.py, so the three
systems compare the same numbers rather than three different random matrices.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import numpy as np

from harness import bench, check, require, seed
from data import spd_int_matrix

require(["numpy.linalg:solve", "numpy.linalg:inv", "numpy.linalg:det",
         "numpy.linalg:eigvals", "numpy.linalg:svd", "numpy.linalg:qr",
         "scipy.linalg:lu"])

seed()
n = 600
a = np.random.random((n, n))
b = np.random.random((n, n))
v = np.random.random(n)
s8 = spd_int_matrix(8)


def r(x, p):
    return int(math.floor(float(x) * 10 ** p + 0.5))


bench("Dot 600x600", lambda: a @ b)
check("Dot 600x600", r((s8 @ s8).sum(), 6))

bench("LinearSolve 600x600", lambda: np.linalg.solve(a, v))
check("LinearSolve 600x600",
      r(np.linalg.solve(s8, np.arange(1.0, 9.0)).sum(), 6))

bench("Inverse 600x600", lambda: np.linalg.inv(a))
check("Inverse 600x600", r(np.linalg.inv(s8).sum(), 6))

bench("Det 600x600", lambda: np.linalg.det(a))
check("Det 600x600", r(np.linalg.det(s8), 4))

bench("Eigenvalues 600x600 (general)", lambda: np.linalg.eigvals(a), reps=1)
check("Eigenvalues 600x600 (general)",
      r(np.sort(np.linalg.eigvals(s8).real).sum(), 4))

m400 = np.random.random((400, 400))
bench("SingularValueDecomposition 400x400", lambda: np.linalg.svd(m400), reps=1)
# Wolfram's SVD part 2 is the diagonal sigma MATRIX; its total is the sum of the
# singular values.
check("SingularValueDecomposition 400x400",
      r(np.linalg.svd(s8, compute_uv=False).sum(), 4))

bench("QRDecomposition 400x400", lambda: np.linalg.qr(m400), reps=1)
# Q and R signs are not canonical across libraries, so the check is |sum(R)|.
check("QRDecomposition 400x400", r(abs(np.linalg.qr(s8)[1].sum()), 4))
