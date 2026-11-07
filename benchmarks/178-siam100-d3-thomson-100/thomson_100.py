"""Experiment 178 -- SIAM 100-Digit Challenge, Appendix D Problem 3 (Python column).

scipy.optimize.minimize (default method: BFGS, finite-difference gradient)
of the numpy-vectorised energy over the same 197 spherical coordinates and
from the same Fibonacci start as thomson_100.m.  mpmath has no minimiser,
so the 50-digit row is a SKIP."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import numpy as np
from scipy.optimize import minimize

require(["scipy.optimize:minimize"])

N = 100
PH1 = np.pi * (3 - np.sqrt(5))
IU = np.triu_indices(N, 1)


def energy(v):
    th = np.concatenate([[0.0], v[:N - 1]])
    ph = np.concatenate([[0.0, PH1], v[N - 1:]])
    X = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], 1)
    d = X[:, None, :] - X[None, :, :]
    return (1.0 / np.sqrt((d ** 2).sum(-1))[IU]).sum()


k = np.arange(N)
v0 = np.concatenate([np.arccos(1 - 2 * k[1:] / N), k[2:] * PH1])
vrun("d3.machine", lambda: minimize(energy, v0).fun, 10)
skip("d3.wp50", "mpmath:minimize")
