"""Experiment 90 -- NMinimize robustness testbed (full corpus), scipy side.

The check carries SOLUTION QUALITY, not the objective value: each case emits
1 if the run reached the published global optimum and 0 if it did not.  A
CHECK-FAIL row therefore means "these two systems disagree about whether they
solved it", which is the robustness gap stated as a measurement instead of a
discarded timing.  See the .m half for the full rationale.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.optimize import differential_evolution, minimize
from harness import bench, check, require

require(["scipy.optimize:differential_evolution", "scipy.optimize:minimize"])

DE = dict(seed=1, maxiter=1500, tol=1e-12, polish=True)
de = lambda f, b: differential_evolution(f, b, **DE)
solved = lambda f, fstar: int(abs(f - fstar) < 0.001)
B = lambda lo, hi, n: [(lo, hi)] * n

def schwefel(v):
    v = np.asarray(v)
    return 418.9828872724339 * v.size - np.sum(v * np.sin(np.sqrt(np.abs(v))))

def griewank(v):
    v = np.asarray(v)
    return (1 + np.sum(v * v) / 4000.0
            - np.prod(np.cos(v / np.sqrt(np.arange(1, v.size + 1)))))

def dropwave(v):
    x, y = v
    return -(1 + np.cos(12 * np.sqrt(x * x + y * y))) / (0.5 * (x * x + y * y) + 2)

def rastrigin(v):
    v = np.asarray(v)
    return 10.0 * v.size + np.sum(v * v - 10.0 * np.cos(2 * np.pi * v))

def styblinski(v):
    v = np.asarray(v)
    return float(0.5 * np.sum(v ** 4 - 16 * v ** 2 + 5 * v))

def bukin6(v):
    x, y = v
    return 100 * np.sqrt(abs(y - 0.01 * x * x)) + 0.01 * abs(x + 10)

def eggholder(v):
    x, y = v
    return (-(y + 47) * np.sin(np.sqrt(abs(y + x / 2 + 47)))
            - x * np.sin(np.sqrt(abs(x - (y + 47)))))

bench("T1 schwefel 5d",        lambda: de(schwefel, B(-500, 500, 5)))
check("T1 schwefel 5d",        solved(de(schwefel, B(-500, 500, 5)).fun, 0.0))

bench("T2 griewank 5d",        lambda: de(griewank, B(-600, 600, 5)))
check("T2 griewank 5d",        solved(de(griewank, B(-600, 600, 5)).fun, 0.0))

bench("T3 drop-wave 2d",       lambda: de(dropwave, B(-5.12, 5.12, 2)))
check("T3 drop-wave 2d",       solved(de(dropwave, B(-5.12, 5.12, 2)).fun, -1.0))

bench("T4 rastrigin 10d",      lambda: de(rastrigin, B(-5.12, 5.12, 10)))
check("T4 rastrigin 10d",      solved(de(rastrigin, B(-5.12, 5.12, 10)).fun, 0.0))

bench("T5 styblinski-tang 5d", lambda: de(styblinski, B(-5, 5, 5)))
check("T5 styblinski-tang 5d", solved(de(styblinski, B(-5, 5, 5)).fun, -195.830828518))

bench("T6 bukin n6",           lambda: de(bukin6, [(-15, -5), (-3, 3)]))
check("T6 bukin n6",           solved(de(bukin6, [(-15, -5), (-3, 3)]).fun, 0.0))

bench("T7 eggholder 2d",       lambda: de(eggholder, B(-512, 512, 2)))
check("T7 eggholder 2d",       solved(de(eggholder, B(-512, 512, 2)).fun, -959.6406627))

# ---- F*: constraint feasibility of the RETURNED point -------------------
q = lambda v: v[0] ** 2 + v[1] ** 2
ineq = [{"type": "ineq", "fun": lambda v: v[0] + v[1] - 2.0}]
eq = [{"type": "eq", "fun": lambda v: v[0] + v[1] - 2.0}]

bench("F1 ineq feasibility",
      lambda: minimize(q, [0.0, 0.0], method="SLSQP", constraints=ineq, tol=1e-12))
_f1 = minimize(q, [0.0, 0.0], method="SLSQP", constraints=ineq, tol=1e-12).x
check("F1 ineq feasibility", int(_f1[0] + _f1[1] >= 2.0))

bench("F2 eq feasibility",
      lambda: minimize(q, [0.0, 0.0], method="SLSQP", constraints=eq, tol=1e-12))
_f2 = minimize(q, [0.0, 0.0], method="SLSQP", constraints=eq, tol=1e-12).x
check("F2 eq feasibility", int(abs(_f2[0] + _f2[1] - 2.0) < 1e-9))
