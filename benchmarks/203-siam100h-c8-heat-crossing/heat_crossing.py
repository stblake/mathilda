"""Experiment 203 -- SIAM 100-Digit Challenge, holdout A for Problem 8 (Python column).

The same ten eigen-series as heat_crossing.m.  scipy has no series summation,
so (as in experiment 173) each series is mpmath.nsum (mpmath.besseljzero for
the disk); the root is scipy.optimize.brentq on a bracket at machine
precision (nsum at dps 15, xtol 1e-15) and mpmath.findroot (secant from the
same start as the .m file) at dps 60.  See heat_crossing.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require, mp
import mpmath
from mpmath import nsum, inf, exp, cos, besseljzero, besselj
from scipy.optimize import brentq

require(["scipy.optimize:brentq", "mpmath:nsum", "mpmath:findroot", "mpmath:besseljzero"])

pi = lambda: mpmath.pi
sg = lambda k: (-1) ** int(k)
w = lambda t: 4 / pi() * nsum(lambda k: sg(k) / (2 * k + 1) * exp(-(2 * k + 1) ** 2 * pi() ** 2 * t / 4), [0, inf])
W = lambda x, t: 4 / pi() * nsum(lambda k: sg(k) / (2 * k + 1) * cos((2 * k + 1) * pi() * x / 2)
                                 * exp(-(2 * k + 1) ** 2 * pi() ** 2 * t / 4), [0, inf])
mean = lambda t: 8 / pi() ** 2 * nsum(lambda k: exp(-(2 * k + 1) ** 2 * pi() ** 2 * t / 4) / (2 * k + 1) ** 2, [0, inf])


def disk(t):
    def term(n):
        j = besseljzero(0, int(n))
        return 2 / (j * besselj(1, j)) * exp(-j * j * t)
    return 1 - nsum(term, [1, inf])


ball = lambda t: 1 + 2 * nsum(lambda n: sg(n) * exp(-n ** 2 * pi() ** 2 * t), [1, inf])
h = mpmath.mpf(1) / 2

# (residual, machine bracket, wp50 start)
C = {
    1: (lambda t: mpmath.mpf(5) / 4 * (1 - w(t) ** 2) - h, (0.05, 1.0), "0.2"),
    2: (lambda t: 1 - w(t) ** 2 - h, (0.05, 1.0), "0.3"),
    3: (lambda t: w(t) - h, (0.05, 1.0), "0.3"),
    4: (lambda t: mean(t) - h, (0.05, 1.0), "0.2"),
    5: (lambda t: 1 - w(t) ** 3 - h, (0.05, 1.0), "0.2"),
    6: (lambda t: 1 - w(t / 4) * w(t) - h, (0.05, 2.0), "0.4"),
    7: (lambda t: disk(t) - h, (0.05, 1.0), "0.15"),
    8: (lambda t: ball(t) - h, (0.02, 1.0), "0.1"),
    9: (lambda t: 1 - W(h, t) * w(t) - h, (0.05, 1.0), "0.2"),
    10: (lambda t: 1 - w(t) ** 2 - mpmath.mpf(1) / 100, (0.01, 0.3), "0.06"),
}


def machine(f, br):
    mp.dps = 15
    return brentq(lambda t: float(f(t)), br[0], br[1], xtol=1e-15, rtol=1e-15)


for i in range(1, 11):
    f, br, s = C[i]
    vrun("hc8.%d machine" % i, lambda f=f, br=br: machine(f, br), 10, 3)
    vrun("hc8.%d wp50" % i, lambda f=f, s=s: mpmath.findroot(f, mpmath.mpf(s)), 48, dps=60)
