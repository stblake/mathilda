"""Experiment 181 -- SIAM 100-Digit Challenge, Appendix D Problem 6 (Python column).

The integral as written over [-1, 1]: scipy.integrate.quad at machine
precision (tight tolerances), mpmath.quad at 60 digits.
See quadruple_exponential.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath
from scipy.integrate import quad

require(["scipy.integrate:quad", "mpmath:quad"])

f = lambda x: math.exp(x + math.sin(math.exp(math.exp(x + 1.0 / 3.0))))
vrun("d6.machine", lambda: quad(f, -1, 1, epsabs=1e-14, epsrel=1e-13, limit=200)[0], 10, 3)

g = lambda x: mpmath.exp(x + mpmath.sin(mpmath.exp(mpmath.exp(x + mpmath.mpf(1) / 3))))
vrun("d6.wp50", lambda: mpmath.quad(g, [-1, 1]), 48, dps=60)
