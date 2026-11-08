"""Experiment 243 -- holdout B for SIAM Appendix D Problem 8 (Python column).

The integrals of rational_rn_integrals.m over all of R^N, as written:
scipy.integrate.quad (N = 1), dblquad (N = 2), tplquad (N = 3) and nquad
(N = 4) on infinite ranges at machine precision; mpmath.quad on
[-inf, inf]^N at 60 digits for N <= 3 (mpmath.quad stops at 3 dimensions,
so the 4-D wp50 row is SKIP)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import math
import mpmath as M
from scipy.integrate import quad, dblquad, tplquad, nquad

require(["scipy.integrate:quad", "scipy.integrate:dblquad", "scipy.integrate:tplquad",
         "scipy.integrate:nquad", "mpmath:quad"])
inf = math.inf
TOL = dict(epsabs=0, epsrel=1e-13)

def mono(N, m):
    """exponent tuples k with |k| <= m"""
    if N == 1:
        return [(i,) for i in range(m + 1)]
    return [(i,) + r for i in range(m + 1) for r in mono(N - 1, m - i)]
def P(ks):
    return lambda *x: sum(math.prod(t ** (2 * e) for t, e in zip(x, k)) for k in ks)
def PM(ks):
    return lambda *x: M.fsum(M.fprod(t ** (2 * e) for t, e in zip(x, k)) for k in ks)

def cmach(N, m):
    p = P(mono(N, m))
    if N == 1:
        I = quad(lambda x: 1 / p(x), -inf, inf, **TOL)[0]
    elif N == 2:
        I = dblquad(lambda y, x: 1 / p(x, y), -inf, inf, -inf, inf, **TOL)[0]
    else:
        I = tplquad(lambda z, y, x: 1 / p(x, y, z), -inf, inf, -inf, inf, -inf, inf, **TOL)[0]
    return math.sqrt(I / (2 * math.pi) ** N)
def cwp(N, m):
    p = PM(mono(N, m))
    I = M.quad(lambda *x: 1 / p(*x), *([[-M.inf, M.inf]] * N))
    return M.sqrt(I / (2 * M.pi) ** N)

for (N, m, k) in [(1, 2, 1), (1, 6, 2), (2, 3, 3), (2, 4, 4), (3, 3, 5), (3, 2, 6)]:
    vrun("bd8.%d machine" % k, lambda: cmach(N, m), 10, 3 if N < 3 else 1, cap=60)
    vrun("bd8.%d wp50" % k, lambda: cwp(N, m), 48, cap=120, dps=60)

vrun("bd8.7 machine", lambda: dblquad(lambda y, x: 1 / (1 + x**4 + y**4), -inf, inf, -inf, inf, **TOL)[0], 10, 3, cap=60)
vrun("bd8.7 wp50", lambda: M.quad(lambda x, y: 1 / (1 + x**4 + y**4), [-M.inf, M.inf], [-M.inf, M.inf]), 48, cap=120, dps=60)
vrun("bd8.8 machine", lambda: tplquad(lambda z, y, x: (1 + x*x + y*y + z*z) ** -3, -inf, inf, -inf, inf, -inf, inf, **TOL)[0], 10, 1, cap=60)
vrun("bd8.8 wp50", lambda: M.quad(lambda x, y, z: (1 + x**2 + y**2 + z**2) ** -3, [-M.inf, M.inf], [-M.inf, M.inf], [-M.inf, M.inf]), 48, cap=120, dps=60)
vrun("bd8.9 machine", lambda: nquad(lambda x, y, z, w: (1 + x*x + y*y + z*z + w*w) ** -3, [[-inf, inf]] * 4,
                                    opts={"epsabs": 0, "epsrel": 1e-13})[0], 10, 1, cap=60)
skip("bd8.9 wp50", "mpmath.quad (at most 3-D)")
vrun("bd8.10 machine", lambda: dblquad(lambda y, x: 1 / (1 + x**4 + x*x*y*y + y**4), -inf, inf, -inf, inf, **TOL)[0], 10, 3, cap=60)
vrun("bd8.10 wp50", lambda: M.quad(lambda x, y: 1 / (1 + x**4 + x**2 * y**2 + y**4), [-M.inf, M.inf], [-M.inf, M.inf]), 48, cap=120, dps=60)
