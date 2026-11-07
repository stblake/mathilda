"""Experiment 175 -- SIAM 100-Digit Challenge, Problem 10 (Python column).

The separated series (4/pi) sum (-1)^k/(2k+1) sech(5(2k+1)pi): scipy has no
series summation, so both rows use mpmath.nsum (dps 15 for the machine row,
dps 60 for wp50).  Python has no PDE solver: c10.machine pde is SKIP.
See brownian_ends.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require, mp
import mpmath

require(["mpmath:nsum"])

p = lambda: 4 / mpmath.pi * mpmath.nsum(
    lambda k: (-1) ** int(k) / (2 * k + 1) * mpmath.sech(5 * (2 * k + 1) * mpmath.pi), [0, mpmath.inf])
vrun("c10.machine", lambda: float(p()), 10, 3, dps=15)
vrun("c10.wp50", p, 48, dps=60)
skip("c10.machine pde", "PDE solver")
