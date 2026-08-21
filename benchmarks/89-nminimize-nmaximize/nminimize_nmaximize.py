"""Experiment 89 -- NMinimize / NMaximize global optimization, timed race.

Competitor side: scipy.optimize.  differential_evolution is the engine-matched
opponent for Mathilda's default DE; minimize(SLSQP) handles the smooth
constrained cases, which is scipy's strongest tool for that shape.

BOUNDS ARE AN ASYMMETRY, NOT A BUG.  scipy's differential_evolution REQUIRES a
box; Mathilda does not, and grows an unbounded coordinate by powers of 10 when
no feasible point is found (nm_driver.c:397-454).  We give scipy each
function's standard published domain, which is the most favourable honest
choice for it -- Mathilda is searching its default +/-10 box with no such hint.
Where the standard domain is tighter than +/-10 this hands scipy an advantage;
that is deliberate and is recorded in the README rather than tuned away.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import math
import numpy as np
from scipy.optimize import differential_evolution, minimize
from harness import bench, check, require

require(["scipy.optimize:differential_evolution", "scipy.optimize:minimize"])

def r(x, p):
    """Match Mathematica's Round[10^p x] on the values this experiment produces."""
    return int(math.floor(x * 10 ** p + 0.5))

SEED = 1
MAXIT = 1500          # mirrors the .m side's MaxIterations -> 1500
DE = dict(seed=SEED, maxiter=MAXIT, tol=1e-12, polish=True)

# ---- objective functions (vector form) ---------------------------------
def rastrigin(v):
    v = np.asarray(v)
    return 10.0 * v.size + np.sum(v * v - 10.0 * np.cos(2 * np.pi * v))

def ackley(v):
    v = np.asarray(v); n = v.size
    return (-20.0 * np.exp(-0.2 * np.sqrt(np.sum(v * v) / n))
            - np.exp(np.sum(np.cos(2 * np.pi * v)) / n) + 20.0 + np.e)

def rosenbrock(v):
    v = np.asarray(v)
    return float(np.sum(100.0 * (v[1:] - v[:-1] ** 2) ** 2 + (1.0 - v[:-1]) ** 2))

def levy(v):
    w = 1.0 + (np.asarray(v) - 1.0) / 4.0
    return float(np.sin(np.pi * w[0]) ** 2
                 + np.sum((w[:-1] - 1) ** 2 * (1 + 10 * np.sin(np.pi * w[:-1] + 1) ** 2))
                 + (w[-1] - 1) ** 2 * (1 + np.sin(2 * np.pi * w[-1]) ** 2))

def sphere(v):
    v = np.asarray(v); return float(np.sum(v * v))

def branin(v):
    x, y = v
    return ((y - 5.1 * x * x / (4 * np.pi ** 2) + 5 * x / np.pi - 6) ** 2
            + 10 * (1 - 1 / (8 * np.pi)) * np.cos(x) + 10)

def camel6(v):
    x, y = v
    return (4 - 2.1 * x * x + x ** 4 / 3) * x * x + x * y + (-4 + 4 * y * y) * y * y

def beale(v):
    x, y = v
    return ((1.5 - x + x * y) ** 2 + (2.25 - x + x * y ** 2) ** 2
            + (2.625 - x + x * y ** 3) ** 2)

def crossintray(v):
    x, y = v
    return -0.0001 * (abs(np.sin(x) * np.sin(y)
                          * np.exp(abs(100 - np.sqrt(x * x + y * y) / np.pi))) + 1) ** 0.1

def styblinski(v):
    v = np.asarray(v)
    return float(0.5 * np.sum(v ** 4 - 16 * v ** 2 + 5 * v))

def dropwave(v):
    x, y = v
    return -(1 + np.cos(12 * np.sqrt(x * x + y * y))) / (0.5 * (x * x + y * y) + 2)

B = lambda lo, hi, n: [(lo, hi)] * n
de = lambda f, b: differential_evolution(f, b, **DE)

# ---- D*: engine-matched DifferentialEvolution --------------------------
bench("D1 rastrigin 2d (DE)",   lambda: de(rastrigin, B(-5.12, 5.12, 2)))
check("D1 rastrigin 2d (DE)",   r(de(rastrigin, B(-5.12, 5.12, 2)).fun, 4))

bench("D2 rastrigin 5d (DE)",   lambda: de(rastrigin, B(-5.12, 5.12, 5)))
check("D2 rastrigin 5d (DE)",   r(de(rastrigin, B(-5.12, 5.12, 5)).fun, 4))

bench("D3 ackley 10d (DE)",     lambda: de(ackley, B(-32.768, 32.768, 10)))
check("D3 ackley 10d (DE)",     r(de(ackley, B(-32.768, 32.768, 10)).fun, 4))

bench("D4 rosenbrock 5d (DE)",  lambda: de(rosenbrock, B(-5.0, 10.0, 5)))
check("D4 rosenbrock 5d (DE)",  r(de(rosenbrock, B(-5.0, 10.0, 5)).fun, 4))

bench("D5 levy 5d (DE)",        lambda: de(levy, B(-10.0, 10.0, 5)))
check("D5 levy 5d (DE)",        r(de(levy, B(-10.0, 10.0, 5)).fun, 4))

bench("D6 sphere 10d (DE)",     lambda: de(sphere, B(-5.12, 5.12, 10)))
check("D6 sphere 10d (DE)",     r(de(sphere, B(-5.12, 5.12, 10)).fun, 4))

# ---- A*: scipy has no "Automatic"; DE is the closest analogue ----------
bench("A1 branin 2d (auto)",       lambda: de(branin, [(-5, 10), (0, 15)]))
check("A1 branin 2d (auto)",       r(de(branin, [(-5, 10), (0, 15)]).fun, 4))

bench("A2 six-hump camel (auto)",  lambda: de(camel6, [(-3, 3), (-2, 2)]))
check("A2 six-hump camel (auto)",  r(de(camel6, [(-3, 3), (-2, 2)]).fun, 4))

bench("A3 beale 2d (auto)",        lambda: de(beale, B(-4.5, 4.5, 2)))
check("A3 beale 2d (auto)",        r(de(beale, B(-4.5, 4.5, 2)).fun, 4))

bench("A4 cross-in-tray (auto)",   lambda: de(crossintray, B(-10, 10, 2)))
check("A4 cross-in-tray (auto)",   r(de(crossintray, B(-10, 10, 2)).fun, 4))

# ---- M*: maximization -- scipy has no maximizer, so negate -------------
# This IS the comparison: Mathilda ships NMaximize, scipy makes you do it by
# hand.  The M3/M4 pair is identical work posed both ways; on the scipy side
# both rows are the same call, so the two timings are a control against which
# Mathilda's wrapper cost (nm_driver.c:555-605) reads directly.
bench("M1 nmaximize styblinski 5d", lambda: de(styblinski, B(-5, 5, 5)))
check("M1 nmaximize styblinski 5d", r(-de(styblinski, B(-5, 5, 5)).fun, 4))

bench("M3 wrapper base nminimize",  lambda: de(rastrigin, B(-5.12, 5.12, 2)))
check("M3 wrapper base nminimize",  r(de(rastrigin, B(-5.12, 5.12, 2)).fun, 4))

# The .m side maximizes -rastrigin, whose maximum is -(min rastrigin) = 0.
# scipy has no maximizer: to maximize g = -rastrigin it minimizes -g =
# +rastrigin and negates the result.  Getting this backwards (minimizing
# -rastrigin, which finds rastrigin's MAXIMUM of ~80.7) is exactly the error
# the value gate caught on the first run of this file -- python reported
# 807066 against 0 from the other two systems.
bench("M4 wrapper same via nmaximize", lambda: de(rastrigin, B(-5.12, 5.12, 2)))
check("M4 wrapper same via nmaximize", r(-de(rastrigin, B(-5.12, 5.12, 2)).fun, 4))

# ---- C*: constrained ---------------------------------------------------
# SLSQP is scipy's strongest smooth-constrained solver.  Rounding is 10^3 to
# match the .m side, which needs it because Mathilda's constrained optimum is
# ~1e-4 infeasible -- see the C* block there and F1/F2 in experiment 90.
q = lambda v: v[0] ** 2 + v[1] ** 2
ineq = [{"type": "ineq", "fun": lambda v: v[0] + v[1] - 2.0}]
eq = [{"type": "eq", "fun": lambda v: v[0] + v[1] - 2.0}]

bench("C1 ineq constrained",
      lambda: minimize(q, [0.0, 0.0], method="SLSQP", constraints=ineq, tol=1e-12))
check("C1 ineq constrained",
      r(minimize(q, [0.0, 0.0], method="SLSQP", constraints=ineq, tol=1e-12).fun, 3))

bench("C2 eq constrained",
      lambda: minimize(q, [0.0, 0.0], method="SLSQP", constraints=eq, tol=1e-12))
check("C2 eq constrained",
      r(minimize(q, [0.0, 0.0], method="SLSQP", constraints=eq, tol=1e-12).fun, 3))

mi = lambda v: (v[0] - 2.4) ** 2 + (v[1] + 1.7) ** 2
mi_kw = dict(seed=SEED, maxiter=MAXIT, tol=1e-12, integrality=[True, True])
bench("C3 mixed integer",
      lambda: differential_evolution(mi, [(-10, 10)] * 2, **mi_kw))
check("C3 mixed integer",
      r(differential_evolution(mi, [(-10, 10)] * 2, **mi_kw).fun, 3))

# ---- I*: indexed-vs-explicit is a Mathilda-only distinction ------------
# scipy sees one vector function in both rows.  Its two timings are therefore
# a control: any divergence between the I1 and I2 ratios is Mathilda's
# indexed-variable dispatch cost and nothing else.
bench("I1 rastrigin 5d explicit vars", lambda: de(rastrigin, B(-5.12, 5.12, 5)))
check("I1 rastrigin 5d explicit vars", r(de(rastrigin, B(-5.12, 5.12, 5)).fun, 4))

bench("I2 rastrigin 5d indexed vars",  lambda: de(rastrigin, B(-5.12, 5.12, 5)))
check("I2 rastrigin 5d indexed vars",  r(de(rastrigin, B(-5.12, 5.12, 5)).fun, 4))
