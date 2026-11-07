"""Experiment 179 -- SIAM 100-Digit Challenge, Appendix D Problem 4 (Python column).

scipy has neither a Moebius function nor an infinite-series summer, so every
row is mpmath: findroot in t = log x (start -34000) on the series as written,
mu(k) from sympy.mobius, li from mpmath.li, the infinite sum by mpmath.nsum --
at mp.dps 15 for the machine row and dps 60 for the 50-digit row; and
findroot on the built-in mpmath.riemannr for the d4.wp50-builtin row.
The answer recorded is x = exp(t*).  See riemann_r_zero.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import mpmath
from sympy import mobius

require(["mpmath:findroot", "mpmath:nsum", "mpmath:li", "sympy:mobius", "mpmath:riemannr"])


def R(t):
    return mpmath.nsum(lambda k: mobius(int(k)) / k * mpmath.li(mpmath.exp(t / k)),
                       [1, mpmath.inf])


vrun("d4.machine", lambda: mpmath.exp(mpmath.findroot(R, -34000)), 10, dps=15)
vrun("d4.wp50", lambda: mpmath.exp(mpmath.findroot(R, -34000)), 48, dps=60)
vrun("d4.wp50-builtin",
     lambda: mpmath.exp(mpmath.findroot(lambda t: mpmath.riemannr(mpmath.exp(t)), -34000)),
     48, dps=60)
