#!/usr/bin/env python3
"""Experiment 53 -- Matrix decompositions (scipy/numpy column).

Same kernels as ``matrix_decompositions.m``, same order and labels.

EXECUTION comparison: numpy/scipy dispatch the same class of LAPACK routine as
Mathilda's src/linalg/ machine path, so a gap is overhead or a missing fast
path, not a missing algorithm.  Baseline is compiled (numpy/scipy over LAPACK).

Timing runs on large random matrices; checks run on the small exactly
reproducible spd_int_matrix(8) or a fixed integer matrix.  Decomposition outputs
are non-unique, so every check is a convention-invariant scalar.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math
import numpy as np
import scipy.linalg as sla

from harness import bench, check, require, bench_if, bench_once, seed
from data import spd_int_matrix

require(["numpy.linalg:svd", "numpy.linalg:qr", "scipy.linalg:lu_factor",
         "numpy.linalg:pinv", "numpy.linalg:matrix_rank", "scipy.linalg:null_space",
         "scipy.linalg:cholesky"])


def r(x, p):
    return int(math.floor(float(x) * 10 ** p + 0.5))


seed()
n = 400
a = np.random.random((n, n))                     # general square, timing
rt = np.random.random((300, 150))                # tall, PseudoInverse timing
w = np.random.random((100, 200))                 # wide, NullSpace timing
m = np.random.random((n, n))
spd = m.T @ m + n * np.eye(n)                     # SPD, Cholesky timing

sm = spd_int_matrix(8)                            # check matrix, well-conditioned
rank_mat = np.array([[i + j for j in range(1, 7)] for i in range(1, 7)],
                    dtype=float)                  # exact rank 2 (outer sum)
rect = np.array([[1.0 / (i + j - 1) for j in range(1, 5)] for i in range(1, 7)])  # 6x4 full col rank

bench("SingularValueDecomposition 400x400", lambda: np.linalg.svd(a), reps=1)
check("SingularValueDecomposition 400x400", r(np.linalg.svd(sm, compute_uv=False).sum(), 6))

bench("QRDecomposition 400x400", lambda: np.linalg.qr(a))
check("QRDecomposition 400x400", r(np.abs(np.diag(np.linalg.qr(sm)[1])).sum(), 6))

bench("LUDecomposition 400x400", lambda: sla.lu_factor(a))
check("LUDecomposition 400x400", r(abs(np.prod(np.diag(sla.lu_factor(sm)[0]))), 4))

bench("PseudoInverse 300x150", lambda: np.linalg.pinv(rt))
check("PseudoInverse 300x150", r(np.linalg.pinv(rect).sum(), 6))

bench("MatrixRank 400x400", lambda: np.linalg.matrix_rank(a))
check("MatrixRank 400x400", int(np.linalg.matrix_rank(rank_mat)))

bench_once("NullSpace 100x200", lambda: sla.null_space(w))
check("NullSpace 100x200", int(sla.null_space(rank_mat).shape[1]))

# CholeskyDecomposition is ABSENT in Mathilda -> the .m side emits SKIP.  scipy
# still measures it, and the single check value passes the gate on its own.
bench_if("CholeskyDecomposition 400x400", "scipy.linalg:cholesky",
         lambda: sla.cholesky(spd), reps=1)
check("CholeskyDecomposition 400x400", r(np.diag(sla.cholesky(sm)).sum(), 6))
