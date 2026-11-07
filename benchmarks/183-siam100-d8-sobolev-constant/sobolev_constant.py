"""Experiment 183 -- SIAM 100-Digit Challenge, Appendix D Problem 8 (Python column).

As sobolev_constant.m: the 10-D integral over R^10 by scipy.integrate.nquad
(scipy's only N-dimensional integrator; integrand through the same
Newton-identity evaluation of Sum_{|k|<=10} y^k) -- expected to hit the cap;
mpmath.quad stops at 3 dimensions, so the 50-digit row is SKIP.  The sanity
row c_3(R^1) is integrate.quad, the reduced row c_2(R^2) integrate.dblquad."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import math
import numpy as np
from scipy.integrate import nquad, quad, dblquad

require(["scipy.integrate:nquad", "scipy.integrate:quad", "scipy.integrate:dblquad"])

def hsum(y):
    p = [sum(v ** j for v in y) for j in range(1, 11)]
    h = [1.0] + [0.0] * 10
    for j in range(1, 11):
        h[j] = sum(p[i - 1] * h[j - i] for i in range(1, j + 1)) / j
    return sum(h)
fS = lambda *x: 1.0 / hsum([t * t for t in x])
inf = math.inf
vrun("d8.machine", lambda: math.sqrt(nquad(fS, [[-inf, inf]] * 10)[0] / (2 * math.pi) ** 10), 3)
skip("d8.wp50", "mpmath.quad (at most 3-D)")
vrun("d8.machine c3R1", lambda: math.sqrt(quad(lambda x: 1 / (1 + x**2 + x**4 + x**6), -inf, inf,
                                                epsabs=0, epsrel=1e-13)[0] / (2 * math.pi)), 10, 3)
vrun("d8.machine c2R2", lambda: math.sqrt(dblquad(lambda y, x: 1 / (1 + x*x + y*y + x**4 + x*x*y*y + y**4),
                                                  -inf, inf, -inf, inf, epsabs=0, epsrel=1e-13)[0] / (2 * math.pi) ** 2), 10, 3)
