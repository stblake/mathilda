"""Experiment 198 -- SIAM 100-digit holdout A, family hc3 (Python column).

The statement of matrix_norms.m (C3's): the 2-norm of the n x m section --
numpy.linalg.norm(A, 2) in doubles, max of mpmath.svd_r (singular values
only) at 60 digits."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import numpy as np
import mpmath as mp

require(["numpy.linalg:norm", "mpmath:svd_r"])

fac = math.factorial
ENT = {1: lambda i, j, F: F(2)**-i * F(3)**-j,
       2: lambda i, j, F: 1/(F(fac(i))*fac(j)),
       3: lambda i, j, F: 1/F(fac(i + j - 1)),
       4: lambda i, j, F: 1/((i + j - 1)*F(2)**(i + j)),
       5: lambda i, j, F: F(1)/(i + j - 1),
       6: lambda i, j, F: F(min(i, j))/max(i, j),
       7: lambda i, j, F: 1/F((i + j - 1)*(i + j - 2)//2 + j)**2,
       8: lambda i, j, F: F(1)/(i + j*j),
       9: lambda i, j, F: F(1)/(1 + abs(i - j))**2,
       10: lambda i, j, F: (mp.cos(i*j) if F is mp.mpf else math.cos(i*j))/(i + j)}
SIZE = {1: (100, 100), 2: (30, 30), 3: (40, 40), 4: (100, 100), 5: (100, 100),
        6: (100, 100), 7: (100, 100), 8: (120, 60), 9: (150, 150), 10: (100, 100)}


def sec(k):
    n, m = SIZE[k]
    return np.array([[ENT[k](i, j, float) for j in range(1, m + 1)] for i in range(1, n + 1)])


def sec_mp(k):
    n, m = SIZE[k]
    return mp.matrix([[ENT[k](i, j, mp.mpf) for j in range(1, m + 1)] for i in range(1, n + 1)])


for k in range(1, 11):
    vrun("hc3.%d machine" % k, lambda: np.linalg.norm(sec(k), 2), 10, 3, cap=60)
    vrun("hc3.%d wp50" % k, lambda: max(mp.svd_r(sec_mp(k), compute_uv=False)), 48, cap=90, dps=60)
