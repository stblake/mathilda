"""Experiment 208 -- SIAM 100-digit holdout A, family d3 (Python column).

scipy.optimize.minimize (default method: BFGS, finite-difference gradient) of
the numpy-vectorised energy over the same coordinates and from the same start
as point_energy_family.m.  mpmath has no minimiser, so every 50-digit row is a
SKIP."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import numpy as np
from scipy.optimize import minimize

require(["scipy.optimize:minimize"])

PH1 = np.pi * (3 - np.sqrt(5))
G = {"c": lambda r2: 1 / np.sqrt(r2), "r2": lambda r2: 1 / r2,
     "log": lambda r2: -np.log(r2) / 2, "dist": lambda r2: -np.sqrt(r2)}
PROBS = [(5, "c"), (12, "c"), (8, "c"), (32, "c"), (20, "r2"), (30, "log"),
         (16, "dist"), (30, "circle"), (50, "c"), (72, "c")]


def make(N, kind):
    IU = np.triu_indices(N, 1)
    k = np.arange(N)
    if kind == "circle":
        def energy(v):
            t = np.concatenate([[0.0], v])
            X = np.stack([np.cos(t), np.sin(t)], 1)
            d = X[:, None, :] - X[None, :, :]
            return (1 / np.sqrt((d ** 2).sum(-1)[IU])).sum()
        return energy, 2 * np.pi * (k[1:] / N) ** 1.5
    g = G[kind]

    def energy(v):
        th = np.concatenate([[0.0], v[:N - 1]])
        ph = np.concatenate([[0.0, PH1], v[N - 1:]])
        X = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], 1)
        d = X[:, None, :] - X[None, :, :]
        return g((d ** 2).sum(-1)[IU]).sum()
    return energy, np.concatenate([np.arccos(1 - 2 * k[1:] / N), k[2:] * PH1])


for i, (N, kind) in enumerate(PROBS, 1):
    energy, v0 = make(N, kind)
    vrun("hd3.%d machine" % i, lambda: minimize(energy, v0).fun, 10, cap=120)
    skip("hd3.%d wp50" % i, "mpmath:minimize")
