"""Experiment 174 -- SIAM 100-Digit Challenge, Problem 9 (Python column).

Machine: the question is a GLOBAL maximum over [0, 5] (the factor
2 + sin 10a makes I(a) multimodal; minimize_scalar(method="bounded") stops at
the local maximum a = 1.3853), so the idiomatic call is the global
scipy.optimize.differential_evolution (seed 0, tol 1e-12, polish) on [0, 5] of
-(2 + sin 10a) * scipy.integrate.quad(x^a sin(a/(2-x)), 0, 2, limit=1000).
wp50: mpmath has no maximiser, so the idiomatic statement is
mpmath.findroot of mpmath.diff(I, a) from a = 0.8, with I(a) = (2 + sin 10a) *
mpmath.quad(x^a sin(a/(2-x)), [0, 2]), at dps 60.  See parametric_maximum.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require, mp
import math
import warnings
warnings.simplefilter("ignore")   # quad's IntegrationWarning on every call
import mpmath
from scipy.integrate import quad
from scipy.optimize import differential_evolution

require(["scipy.optimize:differential_evolution", "scipy.integrate:quad",
         "mpmath:findroot", "mpmath:diff", "mpmath:quad"])


def I(a):
    return (2 + math.sin(10 * a)) * quad(lambda x: x ** a * math.sin(a / (2 - x)), 0, 2, limit=1000)[0]


vrun("c9.machine", lambda: differential_evolution(lambda v: -I(v[0]), [(0, 5)], seed=0,
                                                  tol=1e-12, polish=True).x[0], 10)


def Iw(a):
    return (2 + mpmath.sin(10 * a)) * mpmath.quad(lambda x: x ** a * mpmath.sin(a / (2 - x)), [0, 2])


vrun("c9.wp50", lambda: mpmath.findroot(lambda a: mpmath.diff(Iw, a), mpmath.mpf("0.8")), 48, dps=60)
