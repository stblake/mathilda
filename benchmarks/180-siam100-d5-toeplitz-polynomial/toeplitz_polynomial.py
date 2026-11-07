"""Experiment 180 -- SIAM 100-Digit Challenge, Appendix D Problem 5 (Python column).

scipy.optimize.minimize (default method, BFGS with finite-difference
gradients) of numpy.linalg.norm(p(A), 2) over the 8 lower coefficients from
c = 0, the powers A^j precomputed.  mpmath has no minimiser, so the 50-digit
row is a SKIP.  See toeplitz_polynomial.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import numpy as np
from scipy.optimize import minimize

require(["scipy.optimize:minimize", "numpy.linalg:norm"])

n = 48
A = np.eye(n) - np.eye(n, k=-1) + np.eye(n, k=1) + np.eye(n, k=2) + np.eye(n, k=3)
P = [np.linalg.matrix_power(A, j) for j in range(9)]


def obj(c):
    return np.linalg.norm(P[8] + sum(c[j] * P[j] for j in range(8)), 2)


vrun("d5.machine", lambda: minimize(obj, np.zeros(8)).fun, 10)
skip("d5.wp50", "mpmath:minimize")
