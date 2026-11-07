"""Experiment 171 -- SIAM 100-Digit Challenge, Problem 6 (Python column).

The same statement as biased_flea.m: G(eps) as a 2-D integral over [0,pi]^2
and a root solve of G(eps) = 2 -- scipy.integrate.dblquad (epsabs 0, epsrel
1e-12) inside scipy.optimize.brentq on [0.06, 0.07] at machine precision;
mpmath.quad (2-D tanh-sinh) inside mpmath.findroot (secant from 0.06, 0.07)
at 60 digits."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath
from scipy.integrate import dblquad
from scipy.optimize import brentq

require(["scipy.integrate:dblquad", "scipy.optimize:brentq", "mpmath:quad", "mpmath:findroot"])


def G(e):
    def f(k2, k1):
        c = (math.cos(k1) + math.cos(k2)) / 2
        s = 2 * e * math.sin(k1)
        return (1 - c) / ((1 - c) ** 2 + s * s)
    return dblquad(f, 0, math.pi, 0, math.pi, epsabs=0, epsrel=1e-12)[0] / math.pi ** 2


def Gm(e):
    def f(k1, k2):
        c = (mpmath.cos(k1) + mpmath.cos(k2)) / 2
        s = 2 * e * mpmath.sin(k1)
        return (1 - c) / ((1 - c) ** 2 + s * s)
    return mpmath.quad(f, [0, mpmath.pi], [0, mpmath.pi]) / mpmath.pi ** 2


vrun("c6.machine", lambda: brentq(lambda e: G(e) - 2, 0.06, 0.07, xtol=1e-16, rtol=1e-15), 10, 3)
vrun("c6.wp50", lambda: mpmath.findroot(lambda e: Gm(e) - 2, (mpmath.mpf("0.06"), mpmath.mpf("0.07"))),
     48, dps=60)
