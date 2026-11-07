"""Experiment 210 -- SIAM 100-digit holdout A, family d5 (Python column).

scipy.optimize.minimize (default method: BFGS, finite-difference gradients) of
numpy.linalg.norm(p(A), 2) (or 'fro') over the k free coefficients from c = 0,
the powers A^j precomputed, as experiment 180.  mpmath has no minimiser, so
every 50-digit row is a SKIP.  See matrix_polynomial_family.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import numpy as np
from scipy.optimize import minimize

require(["scipy.optimize:minimize", "numpy.linalg:norm"])


def toep(n, d):
    return sum(v * np.eye(n, k=off) for off, v in d.items())


def grcar(n, s):
    d = {-1: -1.0, 0: 1.0}
    d.update({j: 1.0 for j in range(1, s + 1)})
    return toep(n, d)


PROBS = [
    (np.diag(np.cos(np.arange(9) * np.pi / 8)), 8, "monic", 2),
    (np.diag(np.arange(1.0, 11.0)), 4, "monic", 2),
    (grcar(48, 3), 8, "monic", "fro"),
    (toep(20, {-1: -1.0, 0: 2.0, 1: -1.0}), 5, "monic", 2),
    (grcar(32, 3), 6, "monic", 2),
    (grcar(40, 2), 6, "monic", 2),
    (np.diag(np.arange(1.0, 11.0)), 3, "gmres", 2),
    (toep(12, {0: 1.0, 1: 2.0}), 4, "monic", 2),
    (toep(24, {-1: 0.25, 1: 1.0}), 6, "monic", 2),
    (grcar(64, 3), 10, "monic", 2),
]

for i, (A, k, nz, nm) in enumerate(PROBS, 1):
    P = [np.linalg.matrix_power(A, j) for j in range(k + 1)]
    B0, Bs = (P[k], P[:k]) if nz == "monic" else (P[0], P[1:])

    def obj(c, B0=B0, Bs=Bs, nm=nm):
        return np.linalg.norm(B0 + sum(c[j] * Bs[j] for j in range(len(Bs))), nm)
    vrun("hd5.%d machine" % i, lambda: minimize(obj, np.zeros(k)).fun, 10, cap=60)
    skip("hd5.%d wp50" % i, "mpmath:minimize")
