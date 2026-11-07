"""Experiment 166 -- SIAM 100-Digit Challenge, Problem 1 (Python column).

The integral as written over (0, 1]: scipy.integrate.quad at machine
precision, mpmath.quad at 60 digits.  See twisted_tail.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath
from scipy.integrate import quad

require(["scipy.integrate:quad", "mpmath:quad"])

f = lambda x: math.cos(math.log(x) / x) / x
vrun("c1.machine", lambda: quad(f, 0, 1, limit=1000)[0], 10, 3)

g = lambda x: mpmath.cos(mpmath.log(x) / x) / x
vrun("c1.wp50", lambda: mpmath.quad(g, [0, 1]), 48, dps=60)
