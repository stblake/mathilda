"""Experiment 204 -- SIAM 100-Digit Challenge, holdout A for Problem 9 (Python column).

As in experiment 174.  Machine: the question is a GLOBAL extremum over the
interval, so the idiomatic call is scipy.optimize.differential_evolution
(seed 0, tol 1e-12, polish) of the objective built on scipy.integrate.quad
(limit 1000; plain quad, as a user writes it).
wp50: mpmath has no optimiser, so the statement is mpmath.findroot of
mpmath.diff(F, a) from the .m file's start, F built on mpmath.quad
at dps 60.  Machine rows are capped at 120 s, wp50 rows at 60 s (as
in the .m file).  See parametric_extremum.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require, mp
import math
import warnings
warnings.simplefilter("ignore")
import numpy as np
import mpmath
from mpmath import mpf
from scipy.integrate import quad
from scipy.optimize import differential_evolution

require(["scipy.optimize:differential_evolution", "scipy.integrate:quad",
         "mpmath:findroot", "mpmath:diff", "mpmath:quad"])

Q = lambda f, a, b: quad(f, a, b, limit=1000)[0]
INF = np.inf
# machine objectives (sign: +1 maximise, -1 minimise), interval
FM = {
    1: (lambda a: (2 + math.cos(10 * a)) * Q(lambda x: x ** a * math.sin(a / (2 - x)), 0, 2), 1, (0, 5)),
    2: (lambda a: (2 + math.sin(10 * a)) * Q(lambda x: x ** a * math.sin(a / (1 - x)), 0, 1), 1, (0, 5)),
    3: (lambda a: Q(lambda x: x ** a * math.exp(-x), 0, INF), -1, (0, 2)),
    4: (lambda a: Q(lambda t: math.cos(t - a * math.sin(t)), 0, math.pi) / math.pi, 1, (0, 4)),
    5: (lambda a: (2 + math.sin(3 * a)) * Q(lambda x: x ** -0.5 * math.cos(a * math.log(x)), 0, 1), 1, (0, 4)),
    6: (lambda b: Q(lambda x: math.sin(x * x), 0, b), 1, (0, 3)),
    7: (lambda a: Q(lambda x: x ** a * math.exp(-x) * math.sin(x), 0, INF), 1, (0, 3)),
    8: (lambda a: Q(lambda x: math.exp(-x) * math.cos(a * math.sqrt(x)), 0, INF), -1, (0, 10)),
    9: (lambda a: (2 + math.sin(10 * a)) * Q(lambda x: x ** a * math.cos(a / (2 - x)), 0, 2), 1, (0, 5)),
    10: (lambda a: (2 + math.sin(10 * a)) * Q(lambda x: x ** a * math.sin(a / (3 - x)), 0, 3), 1, (0, 5)),
}

mq = mpmath.quad
ms, mc, me, ml = mpmath.sin, mpmath.cos, mpmath.exp, mpmath.log
FW = {
    1: lambda a: (2 + mc(10 * a)) * mq(lambda x: x ** a * ms(a / (2 - x)), [0, 2]),
    2: lambda a: (2 + ms(10 * a)) * mq(lambda x: x ** a * ms(a / (1 - x)), [0, 1]),
    3: lambda a: mq(lambda x: x ** a * me(-x), [0, mpmath.inf]),
    4: lambda a: mq(lambda t: mc(t - a * ms(t)), [0, mpmath.pi]) / mpmath.pi,
    5: lambda a: (2 + ms(3 * a)) * mq(lambda x: x ** (-mpf(1) / 2) * mc(a * ml(x)), [0, 1]),
    6: lambda b: mq(lambda x: ms(x * x), [0, b]),
    7: lambda a: mq(lambda x: x ** a * me(-x) * ms(x), [0, mpmath.inf]),
    8: lambda a: mq(lambda x: me(-x) * mc(a * mpmath.sqrt(x)), [0, mpmath.inf]),
    9: lambda a: (2 + ms(10 * a)) * mq(lambda x: x ** a * mc(a / (2 - x)), [0, 2]),
    10: lambda a: (2 + ms(10 * a)) * mq(lambda x: x ** a * ms(a / (3 - x)), [0, 3]),
}
START = {1: "0.64", 2: "0.2", 3: "0.46", 4: "1.84", 5: "0.15", 6: "1.77", 7: "1.5", 8: "3", 9: "0.12", 10: "1.42"}

for k in range(1, 11):
    f, sg, (lo, hi) = FM[k]
    vrun("hc9.%d machine" % k, lambda f=f, sg=sg, lo=lo, hi=hi: differential_evolution(
        lambda v: -sg * f(v[0]), [(lo, hi)], seed=0, tol=1e-12, polish=True).x[0], 10, cap=120)
    F = FW[k]
    vrun("hc9.%d wp50" % k, lambda F=F, k=k: mpmath.findroot(lambda a: mpmath.diff(F, a), mpf(START[k])),
         48, cap=60, dps=60)
