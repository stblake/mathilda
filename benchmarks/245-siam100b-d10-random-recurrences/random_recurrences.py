"""Experiment 245 -- holdout B for SIAM Appendix D Problem 10 (Python column).

The Monte-Carlo statements of random_recurrences.m (no library computes a
random-product Lyapunov exponent): machine: n = 10^6 draws of the pair index
from numpy.random.default_rng(1).choice(..., p=weights), ratio recursion
r <- a + b/r in floats, exp(mean(log r)); wp50: exact Python-integer x_n for
n = 10^4, then |x_n|^(1/n) in mpmath at 60 digits."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import numpy as np
import mpmath

require(["numpy.random:default_rng", "mpmath:root"])

C = [([1, 1, 1, 1], [(1, 0), (2, 0), (3, 0), (4, 0)]),
     ([1, 1], [(3, 1), (3, -1)]),
     ([1, 2], [(3, 1), (3, -1)]),
     ([1, 1], [(1, 1), (1, 2)]),
     ([1, 1, 1], [(3, -1), (3, 0), (3, 1)]),
     ([1, 1], [(4, 2), (4, -2)]),
     ([1, 1], [(1, 1), (2, 1)]),
     ([1, 1], [(5, 3), (5, -3)]),
     ([1, 1], [(3, 2), (1, 2)]),
     ([1, 1, 1], [(1, 1), (2, 1), (3, 1)])]
NEED_M = [2, 3, 3, 3, 3, 3, 3, 3, 2, 2]
NEED_W = [1, 2, 2, 2, 2, 2, 2, 2, 1, 1]

def draws(w, n):
    p = np.array(w, dtype=float) / sum(w)
    return np.random.default_rng(1).choice(len(w), n, p=p).tolist()
def lyap_machine(w, cs):
    r, acc = 1.0, 0.0                  # log r_0 = 0
    for i in draws(w, 10 ** 6 - 1):
        a, b = cs[i]
        r = a + b / r
        acc += math.log(r)
    return math.exp(acc / 10 ** 6)
def lyap_wp50(w, cs):
    x0, x1 = 1, 1
    for i in draws(w, 10 ** 4 - 1):
        a, b = cs[i]
        x0, x1 = x1, a * x1 + b * x0
    return mpmath.root(abs(x1), 10 ** 4)

for k, (w, cs) in enumerate(C, 1):
    vrun("bd10.%d machine" % k, lambda: lyap_machine(w, cs), NEED_M[k - 1], 3, cap=60)
    vrun("bd10.%d wp50" % k, lambda: lyap_wp50(w, cs), NEED_W[k - 1], 3, cap=60, dps=60)
