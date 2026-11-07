"""Experiment 200 -- SIAM 100-digit holdout A, family hc5 (Python column).

The statement of discrete_minimax.m (C5's): minimise over the real (and for
k >= 7 imaginary) parts of the coefficients the max over the samples of
|f(x_k) - p(x_k)|, scored against the discrete minimax value.
scipy.optimize.differential_evolution (the global, derivative-free counterpart
of NMinimize; box [-B, B] per coefficient part, B = 2 except k = 3 (3) and
k = 5 (8), fixed seed, tol 1e-12, default polish).  mpmath has no minimiser:
wp50 SKIP."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import numpy as np
from scipy.optimize import differential_evolution
from scipy.special import rgamma, gamma

require(["scipy.optimize:differential_evolution", "scipy.special:rgamma"])

circ = lambda M: np.exp(2j*np.pi*np.arange(M)/M)
PROB = {1: (np.sqrt, np.arange(129)/128, 1, False),
        2: (np.exp, np.cos(np.arange(65)*np.pi/64), 3, False),
        3: (np.abs, -1 + np.arange(129)/64, 4, False),
        4: (rgamma, 2*np.arange(101)/100, 3, False),
        5: (lambda x: np.arctan(5*x), np.cos(np.arange(101)*np.pi/100), 5, False),
        6: (np.log1p, np.arange(65)/64, 2, False),
        7: (rgamma, circ(128), 2, True),
        8: (lambda z: gamma(z + 2), circ(128), 2, True),
        9: (lambda z: np.log(2 + z), circ(64), 2, True),
        10: (rgamma, circ(128), 4, True)}
BOUND = {3: 3, 5: 8}


def solve(k):
    f, x, n, cplx = PROB[k]
    fx = f(x)
    V = np.vander(x, n + 1, increasing=True)
    if cplx:
        obj = lambda c: np.max(np.abs(fx - V @ (c[:n + 1] + 1j*c[n + 1:])))
        nv = 2*(n + 1)
    else:
        obj = lambda c: np.max(np.abs(fx - V @ c))
        nv = n + 1
    B = BOUND.get(k, 2)
    return differential_evolution(obj, [(-B, B)]*nv, tol=1e-12, seed=0, maxiter=20000).fun


for k in range(1, 11):
    vrun("hc5.%d machine" % k, lambda: solve(k), 10, cap=90)
    skip("hc5.%d wp50" % k, "NMinimize")
