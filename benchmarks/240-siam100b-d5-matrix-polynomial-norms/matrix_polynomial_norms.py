"""Experiment 240 -- SIAM 100-Digit Challenge, holdout B for Appendix D
problem 5 (Python column).  scipy.optimize.minimize (default method: BFGS,
finite-difference gradients) of numpy.linalg.norm(p(A), ord) over the k lower
coefficients from c = 0, the powers A^j precomputed -- the same objective and
start as matrix_polynomial_norms.m.  mpmath has no minimiser, so every
50-digit row is a SKIP."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import numpy as np
from scipy.linalg import hilbert
from scipy.optimize import minimize

require(["scipy.optimize:minimize", "numpy.linalg:norm"])


def toep(n, d):
    return sum(v * np.eye(n, k=k) for k, v in d.items())   # key = j - i = numpy offset


def grcar(n):
    return toep(n, {-1: -1, 0: 1, 1: 1, 2: 1, 3: 1})


kac = np.zeros((11, 11))
for j in range(1, 11):
    kac[j - 1, j] = kac[j, j - 1] = np.sqrt(j * (11 - j))


def mcase(lab, A, k, nrm):
    P = [np.linalg.matrix_power(A, j) for j in range(k + 1)]
    obj = lambda c: np.linalg.norm(P[k] + sum(c[j] * P[j] for j in range(k)), nrm)
    vrun(lab + " machine", lambda: minimize(obj, np.zeros(k)).fun, 10, cap=60)
    skip(lab + " wp50", "mpmath:minimize")


mcase("bd5.1", np.diag(np.cos(np.arange(7) * np.pi / 6)), 6, 2)
mcase("bd5.2", grcar(48), 8, "fro")
mcase("bd5.3", grcar(48), 4, 2)
mcase("bd5.4", kac, 4, 2)
mcase("bd5.5", toep(20, {-1: -1, 0: 2, 1: -1}), 4, 2)
mcase("bd5.6", toep(40, {-1: 1, 1: 0.25}), 3, 2)
mcase("bd5.7", grcar(48), 1, 2)
mcase("bd5.8", grcar(24), 6, 2)
mcase("bd5.9", hilbert(8), 3, 2)
mcase("bd5.10", toep(32, {-1: 1, 1: 1, 2: 1}), 3, 2)
