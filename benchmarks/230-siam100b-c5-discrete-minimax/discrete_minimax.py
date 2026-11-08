"""Experiment 230 -- SIAM 100-digit holdout B for Problem 5 (Python column).

The same statement as discrete_minimax.m: minimise over the polynomial's
real coefficients (real problems) or the real and imaginary parts of its
complex coefficients (complex problems) the max over the sample points of
|f - p|, scored against the discrete minimax value E.
scipy.optimize.differential_evolution (the global, derivative-free
counterpart of NMinimize; fixed seed, tol 1e-12, maxiter 20000, default
polish), as in experiment 170.  DE needs a box: each unknown gets
[c - (|c| + 1), c + (|c| + 1)] around the least-squares fit c (numpy lstsq) --
the one piece of glue.  mpmath has no minimiser: wp50 SKIP."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import numpy as np
from scipy.optimize import differential_evolution
from scipy.special import gamma, j0

require(["scipy.optimize:differential_evolution", "scipy.special:gamma", "scipy.special:j0"])

REAL = [(np.exp, (-1, 1), 2, 101), (np.abs, (-1, 1), 4, 201), (np.sqrt, (0, 1), 3, 101),
        (lambda x: 1 / (1 + 25 * x * x), (-1, 1), 6, 201), (gamma, (1, 2), 3, 101),
        (np.arctan, (-1, 1), 5, 201), (lambda x: j0(10 * x), (0, 1), 8, 201)]
CPLX = [(np.exp, 3, 64), (lambda z: np.sqrt(z + 2), 2, 96), (lambda z: gamma(z + 2), 4, 128)]


def setup(k):
    if k < 7:
        f, (a, b), n, M = REAL[k]
        x = a + (b - a) * np.arange(M) / (M - 1)
        V = np.vander(x, n + 1, increasing=True)
        fx = f(x)
        obj = lambda c: np.max(np.abs(fx - V @ c))
        c0 = np.linalg.lstsq(V, fx, rcond=None)[0]
    else:
        f, n, M = CPLX[k - 7]
        z = np.exp(2j * np.pi * np.arange(M) / M)
        V = np.vander(z, n + 1, increasing=True)
        fz = f(z)
        obj = lambda c: np.max(np.abs(fz - V @ (c[:n + 1] + 1j * c[n + 1:])))
        cl = np.linalg.lstsq(V, fz, rcond=None)[0]
        c0 = np.concatenate([cl.real, cl.imag])
    bounds = [(c - abs(c) - 1, c + abs(c) + 1) for c in c0]
    return obj, bounds


def de(k):
    obj, bounds = setup(k)
    return differential_evolution(obj, bounds, tol=1e-12, seed=0, maxiter=20000).fun


for k in range(10):
    vrun("bc5.%d machine" % (k + 1), lambda: de(k), 10, cap=120)
    skip("bc5.%d wp50" % (k + 1), "NMinimize")
