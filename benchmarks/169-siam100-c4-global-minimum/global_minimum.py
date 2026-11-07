"""Experiment 169 -- SIAM 100-Digit Challenge, Problem 4 (Python column).

Global minimiser on the box [-1,1]^2 (see global_minimum.m):
scipy.optimize.differential_evolution (fixed seed, tol 1e-12, default L-BFGS-B
polish) at machine precision.  mpmath has no minimiser, so the wp50 row is that
global search followed by mpmath.findroot on the gradient (mpmath.diff) at 60
digits -- the high-precision polish NMinimize does internally."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath
from scipy.optimize import differential_evolution

require(["scipy.optimize:differential_evolution", "mpmath:findroot"])


def f(v):
    x, y = v
    return (math.exp(math.sin(50 * x)) + math.sin(60 * math.exp(y)) + math.sin(70 * math.sin(x))
            + math.sin(math.sin(80 * y)) - math.sin(10 * (x + y)) + (x * x + y * y) / 4)


def fm(x, y):
    return (mpmath.exp(mpmath.sin(50 * x)) + mpmath.sin(60 * mpmath.exp(y)) + mpmath.sin(70 * mpmath.sin(x))
            + mpmath.sin(mpmath.sin(80 * y)) - mpmath.sin(10 * (x + y)) + (x * x + y * y) / 4)


B = [(-1, 1), (-1, 1)]


def de():
    return differential_evolution(f, B, tol=1e-12, seed=0)


def wp50():
    r = de()
    g = [lambda x, y: mpmath.diff(fm, (x, y), (1, 0)), lambda x, y: mpmath.diff(fm, (x, y), (0, 1))]
    x, y = mpmath.findroot(g, (mpmath.mpf(r.x[0]), mpmath.mpf(r.x[1])))
    return fm(x, y)


vrun("c4.machine", lambda: de().fun, 10, 3)
vrun("c4.wp50", wp50, 48, dps=60)
