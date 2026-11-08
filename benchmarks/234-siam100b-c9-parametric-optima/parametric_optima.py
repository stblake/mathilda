"""Experiment 234 -- SIAM 100-digit HOLDOUT B for Problem 9 (Python column).

Machine: the GLOBAL optimum over the interval, so the idiomatic call is
scipy.optimize.differential_evolution (seed 0, tol 1e-12, polish) of the
objective (negated for an argmax), each objective value one
scipy.integrate.quad call (limit 1000; for the Fourier integrals bc9.3 and
bc9.6 quad's own weight='sin'/'cos' QAWF form, which is how scipy states an
integral to infinity of f(x) sin(w x)).
wp50: mpmath has no optimiser, so the idiomatic statement is mpmath.findroot
of mpmath.diff(F, a) from the .m half's start, F by mpmath.quad on a finite
range and mpmath.quadosc on an infinite oscillatory one (omega = the
frequency; for the Airy integrand bc9.9, whose frequency grows, quadosc's
zeros= form with the asymptotic zeros (3 pi n)^(1/3)), dps 60.
See parametric_optima.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require, mp
import math
import warnings
warnings.simplefilter("ignore")   # quad's IntegrationWarning
import mpmath
from scipy.integrate import quad
from scipy.optimize import differential_evolution

require(["scipy.optimize:differential_evolution", "scipy.integrate:quad",
         "mpmath:findroot", "mpmath:diff", "mpmath:quad", "mpmath:quadosc"])

sin, cos, exp, sqrt, PI = math.sin, math.cos, math.exp, math.sqrt, math.pi
L = dict(limit=1000)
F = {1: lambda a: quad(lambda x: sin(a / x), 0, 1, **L)[0],
     2: lambda a: (2 + sin(10 * a)) * quad(lambda x: x ** a * sin(a / (1 - x)), 0, 1, **L)[0],
     3: lambda a: (2 + cos(7 * a)) * quad(lambda x: x ** -a, 0, math.inf, weight="sin", wvar=1, **L)[0],
     4: lambda a: quad(lambda t: cos(a * t - 10 * sin(t)), 0, PI, **L)[0] / PI,
     5: lambda a: quad(lambda x: exp(-x * x) * sin(a * x), 0, math.inf, **L)[0],
     6: lambda a: quad(lambda x: 1 / (1 + x ** 6), 0, math.inf, weight="cos", wvar=a, **L)[0],
     7: lambda x: quad(lambda t: cos(x * t) / sqrt(1 - t * t), 0, 1, **L)[0],
     8: lambda a: quad(lambda x: sin(a * x * x), 0, 1, **L)[0],
     9: lambda x: quad(lambda t: cos(t ** 3 / 3 + x * t), 0, math.inf, **L)[0] / PI,
     10: lambda b: quad(lambda x: sin(x + b / x), 0, math.inf, **L)[0]}
# (interval, +1 argmax / -1 argmin, wp50 start)
P = {1: ((0, 4), 1, "0.6"), 2: ((0, 5), 1, "0.2"), 3: ((0.1, 1.5), 1, "0.93"), 4: ((0, 5), 1, "2.14"),
     5: ((0, 10), 1, "1.85"), 6: ((0, 6), -1, "4.3"), 7: ((0, 10), -1, "3.83"), 8: ((0, 20), 1, "2.6"),
     9: ((-4, 0), 1, "-1"), 10: ((0.05, 3), 1, "0.2")}

M = mpmath
inf = M.inf
W = {1: lambda a: M.quad(lambda x: M.sin(a / x), [0, 1]),
     2: lambda a: (2 + M.sin(10 * a)) * M.quad(lambda x: x ** a * M.sin(a / (1 - x)), [0, 1]),
     3: lambda a: (2 + M.cos(7 * a)) * M.quadosc(lambda x: x ** -a * M.sin(x), [0, inf], omega=1),
     4: lambda a: M.quad(lambda t: M.cos(a * t - 10 * M.sin(t)), [0, M.pi]) / M.pi,
     5: lambda a: M.quad(lambda x: M.exp(-x * x) * M.sin(a * x), [0, inf]),
     6: lambda a: M.quadosc(lambda x: M.cos(a * x) / (1 + x ** 6), [0, inf], omega=a),
     7: lambda x: M.quad(lambda t: M.cos(x * t) / M.sqrt(1 - t * t), [0, 1]),
     8: lambda a: M.quad(lambda x: M.sin(a * x * x), [0, 1]),
     9: lambda x: M.quadosc(lambda t: M.cos(t ** 3 / 3 + x * t), [0, inf],
                            zeros=lambda n: (3 * M.pi * n) ** (M.mpf(1) / 3)) / M.pi,
     10: lambda b: M.quadosc(lambda x: M.sin(x + b / x), [0, inf], omega=1)}

for k in range(1, 11):
    (lo, hi), sgn, a0 = P[k]
    fk, wk = F[k], W[k]
    vrun("bc9.%d machine" % k, lambda: differential_evolution(lambda v: -sgn * fk(v[0]), [(lo, hi)], seed=0,
                                                             tol=1e-12, polish=True).x[0], 10, cap=85)
    vrun("bc9.%d wp50" % k, lambda: M.findroot(lambda a: M.diff(wk, a), M.mpf(a0)), 48, cap=85, dps=60)
