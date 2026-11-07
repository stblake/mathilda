"""Experiment 173 -- SIAM 100-Digit Challenge, Problem 8 (Python column).

The separated solution u(0,0,t) = (5/4)(1 - w(t)^2) with w the eigen-series
(see heat_plate.m).  scipy has no series summation, so w is mpmath.nsum in
both rows; the root is scipy.optimize.brentq on [0.1, 1] at machine
precision (xtol 1e-15) and mpmath.findroot at dps 60.  scipy has no PDE
solver: c8.machine pde is SKIP.  See heat_plate.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require, mp
import mpmath
from scipy.optimize import brentq

require(["scipy.optimize:brentq", "mpmath:nsum", "mpmath:findroot"])


def w(t):
    return 4 / mpmath.pi * mpmath.nsum(
        lambda k: (-1) ** int(k) / (2 * k + 1) * mpmath.exp(-(2 * k + 1) ** 2 * mpmath.pi ** 2 * t / 4),
        [0, mpmath.inf])


def machine():
    mp.dps = 15
    return brentq(lambda t: float(5 * (1 - w(t) ** 2) / 4 - 1), 0.1, 1.0, xtol=1e-15, rtol=1e-15)


vrun("c8.machine", machine, 10, 3)
vrun("c8.wp50", lambda: mpmath.findroot(lambda t: 5 * (1 - w(t) ** 2) / 4 - 1, mpmath.mpf("0.4")), 48, dps=60)
skip("c8.machine pde", "PDE solver")
