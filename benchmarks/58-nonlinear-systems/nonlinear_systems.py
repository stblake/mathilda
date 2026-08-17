#!/usr/bin/env python3
"""Experiment 58 -- Nonlinear systems (scipy.optimize column).

Same kernels as ``nonlinear_systems.m``, same order and labels.

Baseline is compiled (scipy.optimize.fsolve -> MINPACK hybrd/hybrj).  Each
system has a unique root near the shared start, so both solvers land on the same
solution and the check (a coordinate or the sum of the solution) agrees.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math
import numpy as np
from scipy.optimize import fsolve

from harness import bench, check, require

require(["scipy.optimize:fsolve"])


def r6(x):
    return int(math.floor(float(x) * 1e6 + 0.5))


# ---- small systems ------------------------------------------------------

def poly2(v):
    return [v[0] ** 2 + v[1] ** 2 - 4, v[0] * v[1] - 1]

bench("FindRoot 2x2 polynomial", lambda: fsolve(poly2, [1.5, 0.5]))
check("FindRoot 2x2 polynomial", r6(fsolve(poly2, [1.5, 0.5])[0]))


def transc2(v):
    return [math.cos(v[0]) - v[1], math.sin(v[1]) - v[0] / 2]

bench("FindRoot 2x2 transcendental", lambda: fsolve(transc2, [1.0, 1.0]))
check("FindRoot 2x2 transcendental", r6(fsolve(transc2, [1.0, 1.0]).sum()))


def bf3(v):
    x, y, z = v
    return [3 * x - math.cos(y * z) - 0.5,
            x ** 2 - 81 * (y + 0.1) ** 2 + math.sin(z) + 1.06,
            math.exp(-x * y) + 20 * z + (10 * math.pi - 3) / 3]

bench("FindRoot Burden-Faires 3x3", lambda: fsolve(bf3, [0.1, 0.1, -0.1]))
check("FindRoot Burden-Faires 3x3", r6(fsolve(bf3, [0.1, 0.1, -0.1]).sum()))


# ---- Broyden tridiagonal, N=10 and N=40 ---------------------------------

def broyden_resid(v):
    n = len(v)
    f = np.empty(n)
    for i in range(n):
        left = v[i - 1] if i > 0 else 0.0
        right = v[i + 1] if i < n - 1 else 0.0
        f[i] = (3 - 2 * v[i]) * v[i] - left - 2 * right + 1
    return f


bench("FindRoot Broyden tridiagonal N=10", lambda: fsolve(broyden_resid, -np.ones(10)))
check("FindRoot Broyden tridiagonal N=10", r6(fsolve(broyden_resid, -np.ones(10)).sum()))

bench("FindRoot Broyden tridiagonal N=40", lambda: fsolve(broyden_resid, -np.ones(40)))
check("FindRoot Broyden tridiagonal N=40", r6(fsolve(broyden_resid, -np.ones(40)).sum()))
