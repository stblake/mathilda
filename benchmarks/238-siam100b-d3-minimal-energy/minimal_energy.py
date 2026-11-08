"""Experiment 238 -- SIAM 100-Digit Challenge, holdout B for Appendix D
problem 3 (Python column).  scipy.optimize.minimize (default method: BFGS,
finite-difference gradient) of the numpy-vectorised energy, in the same
coordinates and from the same deterministic starts as minimal_energy.m
(bd3.10 maximises, i.e. minimises the negated sum and reports it negated).
mpmath has no minimiser, so every 50-digit row is a SKIP."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import numpy as np
from scipy.optimize import minimize

require(["scipy.optimize:minimize"])
PH1 = np.pi * (3 - np.sqrt(5))
KERN = {"coul": lambda q: 1 / np.sqrt(q), "log": lambda q: -np.log(q) / 2,
        "r2": lambda q: 1 / q, "dist": lambda q: -np.sqrt(q)}


def sphere(p, N, kd):
    iu = np.triu_indices(N, 1)

    def energy(v):
        th = np.concatenate([[0.0], v[:N - 1]])
        ph = np.concatenate([[0.0, PH1], v[N - 1:]])
        X = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], 1)
        d = X[:, None, :] - X[None, :, :]
        return KERN[kd]((d ** 2).sum(-1)[iu]).sum()
    k = np.arange(N)
    v0 = np.concatenate([np.arccos(1 - 2 * k[1:] / N), k[2:] * PH1])
    sg = -1 if kd == "dist" else 1
    vrun("bd3.%d machine" % p, lambda: sg * minimize(energy, v0).fun, 10, cap=60)
    skip("bd3.%d wp50" % p, "mpmath:minimize")


sphere(1, 5, "coul")
sphere(2, 12, "coul")


def circle(v):
    a = np.concatenate([[0.0], v])
    X = np.stack([np.cos(a), np.sin(a)], 1)
    d = X[:, None, :] - X[None, :, :]
    return (1 / np.sqrt((d ** 2).sum(-1)[np.triu_indices(30, 1)])).sum()


vrun("bd3.3 machine", lambda: minimize(circle, 2 * np.pi * (np.arange(1, 30) / 30) ** 2).fun, 10, cap=60)
skip("bd3.3 wp50", "mpmath:minimize")
sphere(4, 24, "coul")
sphere(5, 32, "coul")
sphere(6, 20, "log")
sphere(7, 20, "r2")


def lj(x):
    X = x.reshape(13, 3)
    d = X[:, None, :] - X[None, :, :]
    q = (d ** 2).sum(-1)[np.triu_indices(13, 1)]
    return (4 * (q ** -6 - q ** -3)).sum()


k = np.arange(12)
t, f = np.arccos(1 - (2 * k + 1) / 12), k * PH1
X0 = np.vstack([[0, 0, 0], 1.1 * np.stack([np.sin(t) * np.cos(f), np.sin(t) * np.sin(f), np.cos(t)], 1)])
vrun("bd3.8 machine", lambda: minimize(lj, X0.ravel()).fun, 10, cap=60)
skip("bd3.8 wp50", "mpmath:minimize")
sphere(9, 60, "coul")
sphere(10, 30, "dist")
