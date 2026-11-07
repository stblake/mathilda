"""Experiment 199 -- SIAM 100-digit holdout A, family hc4 (Python column).

Global minimiser on the box (see global_minima.m), C4's statement:
scipy.optimize.differential_evolution (fixed seed, tol 1e-12, default L-BFGS-B
polish) at machine precision.  mpmath has no minimiser, so the wp50 row is that
global search followed by mpmath.findroot on the gradient (mpmath.diff) at 60
digits -- the high-precision polish NMinimize does internally (for hc4.8, whose
minimiser lies on the box edge, the gradient is not zero there: a finding for
this formulation, not adjusted)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath as mp
from scipy.optimize import differential_evolution

require(["scipy.optimize:differential_evolution", "mpmath:findroot"])


def objs(L):
    s, c, e, sq, ab, pi = L.sin, L.cos, L.exp, L.sqrt, L.fabs, L.pi
    one = mp.mpf(1) if L is mp else 1.0
    return {
        1: lambda x, y: (4 - 21*one/10*x**2 + x**4/3)*x**2 + x*y + (-4 + 4*y**2)*y**2,
        2: lambda x, y: (y - 51*one/10/(4*pi**2)*x**2 + 5/pi*x - 6)**2 + 10*(1 - 1/(8*pi))*c(x) + 10,
        3: lambda x, y, z: sum((u**4 - 16*u**2 + 5*u)/2 for u in (x, y, z)),
        4: lambda x, y: sum(i*c((i + 1)*x + i) for i in range(1, 6))*sum(i*c((i + 1)*y + i) for i in range(1, 6)),
        5: lambda x, y: -ab(s(x)*c(y)*e(ab(1 - sq(x*x + y*y)/pi))),
        6: lambda x, y: -(ab(s(x)*s(y)*e(ab(100 - sq(x*x + y*y)/pi))) + 1)**(one/10)/10000,
        7: lambda x, y: e(s(30*x)) + s(40*e(y)) + s(50*s(x)) + s(s(60*y)) - s(10*(x + y)) + (x*x + y*y)/4,
        8: lambda x, y: -(y + 47)*s(sq(ab(x/2 + y + 47))) - x*s(sq(ab(x - (y + 47)))),
        9: lambda x, y: (x*x + y - 11)**2 + (x + y*y - 7)**2 + x + y,
        10: lambda x, y: e(s(70*x)) + s(80*e(y)) + s(90*s(x)) + s(s(100*y)) - s(12*(x + y)) + (x*x + y*y)/4}


FM, FP = objs(math), objs(mp)
BOX = {1: [(-3, 3), (-2, 2)], 2: [(-5, 10), (0, 15)], 3: [(-5, 5)]*3, 4: [(-10, 10)]*2,
       5: [(-10, 10)]*2, 6: [(-10, 10)]*2, 7: [(-1, 1)]*2, 8: [(-512, 512)]*2,
       9: [(-5, 5)]*2, 10: [(-1, 1)]*2}


def de(k):
    return differential_evolution(lambda v: FM[k](*v), BOX[k], tol=1e-12, seed=0)


def wp50(k):
    r = de(k)
    n = len(BOX[k])
    g = [(lambda i: (lambda *v: mp.diff(FP[k], v, tuple(1 if j == i else 0 for j in range(n)))))(i)
         for i in range(n)]
    v = mp.findroot(g, tuple(mp.mpf(t) for t in r.x))
    return FP[k](*(v if n > 1 else [v]))


for k in range(1, 11):
    vrun("hc4.%d machine" % k, lambda: de(k).fun, 10, cap=60)
    vrun("hc4.%d wp50" % k, lambda: wp50(k), 48, cap=90, dps=60)
