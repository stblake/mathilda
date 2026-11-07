"""Experiment 228 -- SIAM 100-digit holdout B for Problem 3 (Python column).

The same statement as matrix_norms.m: the 2-norm of each section --
numpy.linalg.norm(A, 2) in doubles, mpmath.svd_r / svd_c (singular values
only) at 60 digits."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp

require(["numpy.linalg:norm", "mpmath:svd_r", "mpmath:svd_c"])


def entry(k, i, j, n, num):
    """a_ij of matrix k (1-based i, j) in the number type num (float / mpf)."""
    if k == 1: return num(1) / (i + j - 1)
    if k == 2: return num(1) / mp.factorial(i + j - 1) if num is mp.mpf else 1.0 / float(mp.factorial(i + j - 1))
    if k == 3: return num(1) / (1 + (i - j) ** 2)
    if k == 4: return num(1) / i if j <= i else num(0)
    if k == 5: return num(1) / 2 ** (i - j) if i >= j else num(-1) ** (j - i) / num(3) ** (j - i)
    if k == 6: return num(1) / ((i + j - 1) * (i + j - 2) // 2 + j)
    if k == 7:
        r = (i * j) % 7                                  # exact phase reduction
        if num is mp.mpf: return mp.expjpi(mp.mpf(2 * r) / 7) / (i + j)
        return np.exp(2j * np.pi * r / 7) / (i + j)
    if k == 8: return num(n + 1 - max(i, j)) if j >= i - 1 else num(0)
    if k == 9: return num(1) / (num(2) ** i + num(3) ** j)
    if k == 10: return 1 / mp.sqrt(i + j - 1) if num is mp.mpf else 1 / np.sqrt(i + j - 1)


def build(k, size, num):
    m, n = size if isinstance(size, tuple) else (size, size)
    rows = [[entry(k, i, j, n, num) for j in range(1, n + 1)] for i in range(1, m + 1)]
    return np.array(rows) if num is float else mp.matrix(rows)


def norm_mach(k, size):
    return np.linalg.norm(build(k, size, float), 2)


def norm_mp(k, size):
    A = build(k, size, mp.mpf)
    sv = mp.svd_c(A, compute_uv=False) if k == 7 else mp.svd_r(A, compute_uv=False)
    return max(sv)


SIZES = {1: (500, 60), 2: (40, 40), 3: (400, 60), 4: (500, 80), 5: (400, 60),
         6: ((500, 150), (100, 30)), 7: (300, 50), 8: (300, 60), 9: (150, 150), 10: (1000, 120)}

for k in range(1, 11):
    nm, nw = SIZES[k]
    vrun("bc3.%d machine" % k, lambda: norm_mach(k, nm), 10, 3, cap=60)
    vrun("bc3.%d wp50" % k, lambda: norm_mp(k, nw), 48, cap=120, dps=60)
