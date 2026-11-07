"""Experiment 168 -- SIAM 100-Digit Challenge, Problem 3 (Python column).

The same statement as operator_norm.m: the 2-norm of the n x n section A_n
-- numpy.linalg.norm(A, 2) in doubles, mpmath.svd_r (singular values only)
at 60 digits.  n = 4096 for the full problem, 500 / 100 for the reduced
rows."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import numpy as np
import mpmath

require(["numpy.linalg:norm", "mpmath:svd_r"])


def sec(n):
    i = np.arange(1, n + 1)[:, None]
    j = np.arange(1, n + 1)[None, :]
    return 1.0 / ((i + j - 1) * (i + j - 2) / 2 + j)


def sec_mp(n):
    return mpmath.matrix([[mpmath.mpf(1) / ((i + j - 1) * (i + j - 2) // 2 + j)
                           for j in range(1, n + 1)] for i in range(1, n + 1)])


def norm_mp(n):
    return max(mpmath.svd_r(sec_mp(n), compute_uv=False))


vrun("c3.machine", lambda: np.linalg.norm(sec(4096), 2), 10)
vrun("c3.wp50", lambda: norm_mp(4096), 38, cap=60, dps=60)
vrun("c3.machine-n500", lambda: np.linalg.norm(sec(500), 2), 10, 3)
vrun("c3.wp50-n100", lambda: norm_mp(100), 48, dps=60)
