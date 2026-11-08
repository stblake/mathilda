"""Experiment 229 -- SIAM 100-digit holdout B for Problem 4 (Python column).

Global minimiser on each box (see global_minima.m):
scipy.optimize.differential_evolution (fixed seed, tol 1e-12, default L-BFGS-B
polish) at machine precision.  mpmath has no minimiser, so the wp50 row is that
global search followed by mpmath.findroot on the gradient (mpmath.diff) at 60
digits -- the high-precision polish NMinimize does internally, exactly as in
experiment 169.  (bc4.6 has its minimum on the boundary, where the gradient is
not zero: that polish then goes elsewhere -- a limit of the formulation, kept.)"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath as mp
from scipy.optimize import differential_evolution

require(["scipy.optimize:differential_evolution", "mpmath:findroot"])


def make(lib):
    """The ten objectives over the number library lib (math or mpmath)."""
    sin, cos, exp, sqrt = lib.sin, lib.cos, lib.exp, lib.sqrt
    pi, E = (math.pi, math.e) if lib is math else (mp.pi, mp.e)
    c = (lambda a, b: a / b) if lib is math else (lambda a, b: mp.mpf(a) / b)
    g5 = lambda t: sum(i * cos((i + 1) * t + i) for i in range(1, 6))
    return [
        lambda x, y: sin(10 * x) + cos(11 * y) + (x * x + y * y) / 4,
        lambda x, y: (x * x + y - 11) ** 2 + (x + y * y - 7) ** 2 + 3 * x,
        lambda x, y, z: (exp(sin(30 * x)) + sin(35 * exp(y)) + sin(sin(40 * z))
                         - sin(10 * (x + y + z)) + (x * x + y * y + z * z) / 4),
        lambda x, y: (-20 * exp(-sqrt(((x - c(3, 10)) ** 2 + (y + c(2, 5)) ** 2) / 2) / 5)
                      - exp((cos(2 * pi * x) + cos(2 * pi * y)) / 2) + E + 20),
        lambda x, y: g5(x) * g5(y),
        lambda x, y: (-(y + 47) * sin(sqrt(abs(x / 2 + y + 47)))
                      - x * sin(sqrt(abs(x - (y + 47))))),
        lambda x: cos(x) + cos(sqrt(2) * x) + cos(sqrt(3) * x),
        lambda x, y: ((sin(25 * x) + sin(25 * y)) / 2 + (x * x + y * y) / 10
                      - 2 * exp(-300 * ((x - c(3, 5)) ** 2 + (y + c(7, 20)) ** 2))),
        lambda x, y: 100 * (y - x * x) ** 2 + (1 - x) ** 2 + c(3, 10) * sin(15 * x) * sin(15 * y),
        lambda x, y, z, w: (sum(v * v - 10 * cos(2 * pi * v) for v in (x, y, z, w))
                            + (x * y + y * z + z * w) / 2 + c(3, 10) * x),
    ]


FM, FP = make(math), make(mp)
BOX = [[(-1, 1)] * 2, [(-5, 5)] * 2, [(-1, 1)] * 3, [(-5, 5)] * 2, [(-10, 10)] * 2,
       [(-512, 512)] * 2, [(0, 100)], [(-1, 1)] * 2, [(-2, 2)] * 2, [(-5, 5)] * 4]


def de(k):
    return differential_evolution(lambda v: FM[k](*v), BOX[k], tol=1e-12, seed=0)


def wp50(k):
    r = de(k)
    d = len(BOX[k])
    f = FP[k]
    grad = [(lambda i: (lambda *v: mp.diff(f, v, tuple(int(j == i) for j in range(d)))))(i)
            for i in range(d)]
    x0 = [mp.mpf(float(t)) for t in r.x]
    if d == 1:
        v = [mp.findroot(grad[0], x0[0])]
    else:
        v = list(mp.findroot(grad, x0))
    return f(*v)


for k in range(10):
    vrun("bc4.%d machine" % (k + 1), lambda: de(k).fun, 10, 3, cap=60)
    vrun("bc4.%d wp50" % (k + 1), lambda: wp50(k), 48, cap=120, dps=60)
