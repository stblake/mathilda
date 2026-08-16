#!/usr/bin/env python3
"""Experiment 79 -- SHGO global optimization (scipy column).

Same 8 problems as ``shgo.m``, same labels, same order, solved with
scipy.optimize.shgo using the SAME sampling_method and sample count (n) on both
sides, so the race is SHGO-vs-shgo (both build a sampled complex, find the
minimizer pool, and multi-start a local minimizer from it). These are the
standard bounded multimodal benchmarks -- Himmelblau, Branin, six-hump camel,
Cross-in-tray, Shubert, Rastrigin, Ackley, Styblinski-Tang -- where the pool +
homology construction is designed to shine. The check is the objective at the
global optimum rounded to 10^6, which both systems reach and polish to
agreement.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import shgo

from harness import bench, check, require

require(["scipy.optimize:shgo"])


def _solve(fun, bounds, n, method, iters=1):
    return shgo(fun, bounds, n=n, iters=iters, sampling_method=method)


# ---- 2-D benchmarks (simplicial) ----
def _himmelblau(v):
    x, y = v
    return (x * x + y - 11.0) ** 2 + (x + y * y - 7.0) ** 2

def _branin(v):
    x, y = v
    return ((y - 5.1 / (4 * np.pi ** 2) * x * x + 5 / np.pi * x - 6.0) ** 2
            + 10 * (1 - 1 / (8 * np.pi)) * np.cos(x) + 10)

def _camel(v):
    x, y = v
    return (4 - 2.1 * x * x + x ** 4 / 3) * x * x + x * y + (-4 + 4 * y * y) * y * y

def _cross_in_tray(v):
    x, y = v
    return -0.0001 * (abs(np.sin(x) * np.sin(y)
                          * np.exp(abs(100 - np.sqrt(x * x + y * y) / np.pi))) + 1) ** 0.1

def _shubert(v):
    x, y = v
    sx = sum(j * np.cos((j + 1) * x + j) for j in range(1, 6))
    sy = sum(j * np.cos((j + 1) * y + j) for j in range(1, 6))
    return sx * sy


# ---- n-D benchmarks ----
def _rastrigin(v):
    v = np.asarray(v)
    return 10.0 * v.size + np.sum(v * v - 10.0 * np.cos(2 * np.pi * v))

def _ackley(v):
    v = np.asarray(v)
    n = v.size
    return (-20.0 * np.exp(-0.2 * np.sqrt(np.sum(v * v) / n))
            - np.exp(np.sum(np.cos(2 * np.pi * v)) / n) + 20.0 + np.e)

def _styblinski(v):
    v = np.asarray(v)
    return float(np.sum(v ** 4 - 16 * v * v + 5 * v) / 2.0)


_c01 = lambda: _solve(_himmelblau,   [(-5, 5), (-5, 5)],          100, "simplicial")
_c02 = lambda: _solve(_branin,       [(-5, 10), (0, 15)],         100, "simplicial")
_c03 = lambda: _solve(_camel,        [(-3, 3), (-2, 2)],          200, "simplicial")
_c04 = lambda: _solve(_cross_in_tray,[(-10, 10), (-10, 10)],      500, "simplicial")
_c05 = lambda: _solve(_shubert,      [(-10, 10), (-10, 10)],      600, "simplicial")
_c06 = lambda: _solve(_rastrigin,    [(-5.12, 5.12)] * 3,         200, "simplicial")
_c07 = lambda: _solve(_ackley,       [(-5, 5)] * 4,               200, "simplicial")
_c08 = lambda: _solve(_styblinski,   [(-5, 5)] * 4,               256, "sobol")

bench("01 Himmelblau simplicial", _c01)
check("01 Himmelblau simplicial", int(np.floor(1e6 * _c01().fun + 0.5)))

bench("02 Branin simplicial", _c02)
check("02 Branin simplicial", int(np.floor(1e6 * _c02().fun + 0.5)))

bench("03 Six-hump camel simplicial", _c03)
check("03 Six-hump camel simplicial", int(np.floor(1e6 * _c03().fun + 0.5)))

bench("04 Cross-in-tray simplicial", _c04)
check("04 Cross-in-tray simplicial", int(np.floor(1e6 * _c04().fun + 0.5)))

bench("05 Shubert simplicial", _c05)
check("05 Shubert simplicial", int(np.floor(1e6 * _c05().fun + 0.5)))

bench("06 Rastrigin 3D simplicial", _c06)
check("06 Rastrigin 3D simplicial", int(np.floor(1e6 * _c06().fun + 0.5)))

bench("07 Ackley 4D simplicial", _c07)
check("07 Ackley 4D simplicial", int(np.floor(1e6 * _c07().fun + 0.5)))

bench("08 Styblinski-Tang 4D sobol", _c08)
check("08 Styblinski-Tang 4D sobol", int(np.floor(1e6 * _c08().fun + 0.5)))
