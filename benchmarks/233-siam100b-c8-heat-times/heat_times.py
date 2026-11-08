"""Experiment 233 -- SIAM 100-digit HOLDOUT B for Problem 8 (Python column).

The same ten statements as heat_times.m: each solution as its separated
series summed by mpmath.nsum (scipy has no series summation -- as in
experiment 173), the root by scipy.optimize.brentq on a bracket around it at
machine precision (nsum at dps 15) and by mpmath.findroot (secant from the
same starting value as the .m half) at dps 60.  The Bessel zeros of bc8.4 are
mpmath.besseljzero.  See heat_times.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require, mp
import mpmath
from mpmath import pi, exp, nsum, inf, sin, erfc, sqrt, besseljzero, besselj
from scipy.optimize import brentq

require(["scipy.optimize:brentq", "mpmath:nsum", "mpmath:findroot", "mpmath:besseljzero"])

sg = lambda k: -1 if int(k) % 2 else 1
w = lambda t: 4 / pi * nsum(lambda k: sg(k) / (2 * k + 1) * exp(-(2 * k + 1) ** 2 * pi ** 2 * t / 4), [0, inf])
w1 = lambda t: 4 / pi * nsum(lambda k: sg(k) / (2 * k + 1) * exp(-(2 * k + 1) ** 2 * pi ** 2 * t), [0, inf])
u2 = lambda t: mpmath.mpf(1) / 2 - 2 / pi * nsum(lambda k: sg(k) / (2 * k + 1) * exp(-(2 * k + 1) ** 2 * pi ** 2 * t), [0, inf])


def u4(t):
    j = lambda k: besseljzero(0, int(k))
    return 1 - 2 * nsum(lambda k: exp(-j(k) ** 2 * t) / (j(k) * besselj(1, j(k))), [1, inf])


u5 = lambda t: 1 + 2 * nsum(lambda n: sg(n) * exp(-n ** 2 * pi ** 2 * t), [1, inf])
wb = lambda t: 8 / pi ** 2 * nsum(lambda k: exp(-(2 * k + 1) ** 2 * pi ** 2 * t / 4) / (2 * k + 1) ** 2, [0, inf])
u8 = lambda t: 2 * nsum(lambda n: sg(n) * erfc((2 * n + 1) / (4 * sqrt(t))), [0, inf])
v9 = lambda t: 4 / pi * nsum(lambda k: sin((2 * k + 1) * pi / 4) / (2 * k + 1) * exp(-(2 * k + 1) ** 2 * pi ** 2 * t), [0, inf])

F = mpmath.mpf
# label: (residual g(t) whose root is asked, bracket for brentq, start for findroot)
P = [("bc8.1", lambda t: F(5) / 4 * (1 - w(t) ** 2) - F(3) / 10, (0.1, 0.2), "0.15"),
     ("bc8.2", lambda t: u2(t) - F(1) / 4, (0.05, 0.15), "0.1"),
     ("bc8.3", lambda t: 1 - w(t) ** 3 - F(9) / 10, (0.3, 0.5), "0.4"),
     ("bc8.4", lambda t: u4(t) - F(1) / 2, (0.1, 0.3), "0.2"),
     ("bc8.5", lambda t: u5(t) - F(1) / 2, (0.1, 0.2), "0.14"),
     ("bc8.6", lambda t: F(5) / 4 * (1 - wb(t) ** 2) - 1, (0.15, 0.35), "0.24"),
     ("bc8.7", lambda t: 1 - w(t) * w1(t) - F(1) / 2, (0.05, 0.15), "0.09"),
     ("bc8.8", lambda t: u8(t) - F(1) / 10, (0.02, 0.05), "0.03"),
     ("bc8.9", lambda t: 1 - v9(t) ** 2 - F(1) / 2, (0.02, 0.04), "0.03"),
     ("bc8.10", lambda t: F(5) / 4 * (1 - w(t) ** 2) - F(10) ** -6, (0.012, 0.025), "0.02")]


def machine(g, br):
    def run():
        mp.dps = 15
        return brentq(lambda t: float(g(t)), br[0], br[1], xtol=1e-16, rtol=1e-15)
    return run


for lab, g, br, t0 in P:
    vrun(lab + " machine", machine(g, br), 10, cap=80)
    vrun(lab + " wp50", lambda: mpmath.findroot(g, F(t0)), 48, cap=80, dps=60)
