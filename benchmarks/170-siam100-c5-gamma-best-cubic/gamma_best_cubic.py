"""Experiment 170 -- SIAM 100-Digit Challenge, Problem 5 (Python column).

The same statement as gamma_best_cubic.m: minimise over the 8 real and
imaginary parts of c0..c3 the max over the 128th roots of unity z_k of
|1/Gamma(z_k) - (c0 + c1 z + c2 z^2 + c3 z^3)|, scored against the discrete
minimax value E_128.  scipy.optimize.differential_evolution (the global,
derivative-free counterpart of NMinimize; box [-2,2]^8, fixed seed,
tol 1e-12, default polish) with scipy.special.rgamma.  mpmath has no
minimiser: wp50 SKIP."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import numpy as np
from scipy.optimize import differential_evolution
from scipy.special import rgamma

require(["scipy.optimize:differential_evolution", "scipy.special:rgamma"])

M = 128
zk = np.exp(2j * np.pi * np.arange(M) / M)
fk = rgamma(zk)
V = np.vander(zk, 4, increasing=True)


def obj(x):
    return np.max(np.abs(fk - V @ (x[:4] + 1j * x[4:])))


vrun("c5.machine", lambda: differential_evolution(obj, [(-2, 2)] * 8, tol=1e-12, seed=0,
                                                  maxiter=20000).fun, 10)
skip("c5.wp50", "NMinimize")
