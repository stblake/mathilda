"""Experiment 185 -- SIAM 100-Digit Challenge, Appendix D Problem 10 (Python column).

The Monte-Carlo statement of random_recurrence.m (no library computes a
random-product Lyapunov exponent):  machine: n = 10^6 signs from
numpy.random.default_rng(1), ratio recursion r <- 2 +- 1/r in floats,
exp(mean(log r)); wp50: exact Python-integer x_n for n = 10^4, then
|x_n|^(1/n) in mpmath at 60 digits."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import numpy as np
import mpmath

require(["numpy.random:default_rng", "mpmath:root"])

def lyap_machine():
    s = np.random.default_rng(1).choice([-1.0, 1.0], 10 ** 6 - 1).tolist()
    r, acc = 1.0, 0.0                  # log r_0 = 0
    for e in s:
        r = 2.0 + e / r
        acc += math.log(r)
    return math.exp(acc / 10 ** 6)
vrun("d10.machine", lyap_machine, 2, 3)

def lyap_wp50():
    s = np.random.default_rng(1).choice([-1, 1], 10 ** 4 - 1).tolist()
    a, b = 1, 1
    for e in s:
        a, b = b, 2 * b + e * a
    return mpmath.root(abs(b), 10 ** 4)
vrun("d10.wp50", lyap_wp50, 1, 3, dps=60)
