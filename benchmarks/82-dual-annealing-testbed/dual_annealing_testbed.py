#!/usr/bin/env python3
"""Experiment 82 -- Dual Annealing across the NMinimize test corpus (scipy column).

Same 9 problems as ``dual_annealing_testbed.m``, same labels, same order, solved
with scipy.optimize.dual_annealing on the identical bounded box, scipy's default
parameters and a matched seed=1 on both sides. This is the cleanly-comparable
subset of the NMinimize test corpus -- the continuous, box-bounded multimodal
stress functions where both implementations reach the same global. The rest of
the corpus is analysed in README.md: scipy.optimize.dual_annealing is BOX-ONLY
(no constraints argument), so the constrained / equality / disjunctive /
combinatorial / mixed-integer cases are not dual-annealing races, and a handful of
sharp-needle functions are missed by both.

The check is the objective at the global optimum, rounded per problem so both
implementations, reaching the same basin, agree.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import dual_annealing

from harness import bench, check, require

require(["scipy.optimize:dual_annealing"])

SEED = 1


def _solve(fun, bounds):
    return dual_annealing(fun, bounds, seed=SEED)


def _val(res, scale):
    return int(np.floor(scale * res.fun + 0.5))


def _chained(v):
    x = v[0]
    return np.sin(2 * x) + np.cos(x)

def _schwefel2(v):
    x, y = v
    return 837.9658 - x * np.sin(np.sqrt(abs(x))) - y * np.sin(np.sqrt(abs(y)))

def _schaffer2(v):
    x, y = v
    return 0.5 + (np.sin(x * x - y * y) ** 2 - 0.5) / (1 + 0.001 * (x * x + y * y)) ** 2

def _rugged(v):
    x, y = v
    return x * x + y * y + 10 * np.sin(3 * x) ** 2 + 10 * np.sin(3 * y) ** 2

def _ackley(v):
    v = np.asarray(v)
    n = v.size
    return (-20 * np.exp(-0.2 * np.sqrt(np.sum(v * v) / n))
            - np.exp(np.sum(np.cos(2 * np.pi * v)) / n) + 20 + np.e)

def _rastrigin(v):
    v = np.asarray(v)
    return 10 * v.size + np.sum(v * v - 10 * np.cos(2 * np.pi * v))

def _rosenbrock(v):
    v = np.asarray(v)
    return float(np.sum(100 * (v[1:] - v[:-1] ** 2) ** 2 + (1 - v[:-1]) ** 2))

def _schwefel(v):
    v = np.asarray(v)
    return 418.9829 * v.size - np.sum(v * np.sin(np.sqrt(np.abs(v))))


_c01 = lambda: _solve(_chained,    [(-2, 3)])
_c02 = lambda: _solve(_schwefel2,  [(-500, 500)] * 2)
_c03 = lambda: _solve(_schaffer2,  [(-100, 100)] * 2)
_c04 = lambda: _solve(_rugged,     [(-5, 5)] * 2)
_c05 = lambda: _solve(_ackley,     [(-32, 32)] * 3)
_c06 = lambda: _solve(_rastrigin,  [(-5.12, 5.12)] * 5)
_c07 = lambda: _solve(_rastrigin,  [(-5.12, 5.12)] * 8)
_c08 = lambda: _solve(_rosenbrock, [(-5, 5)] * 10)
_c09 = lambda: _solve(_schwefel,   [(-500, 500)] * 10)

bench("01 Chained trig 1D", _c01)
check("01 Chained trig 1D", _val(_c01(), 1e4))

bench("02 Schwefel 2D", _c02)
check("02 Schwefel 2D", _val(_c02(), 1e3))

bench("03 Schaffer N2 2D", _c03)
check("03 Schaffer N2 2D", _val(_c03(), 1e4))

bench("04 Rugged sine 2D", _c04)
check("04 Rugged sine 2D", _val(_c04(), 1e4))

bench("05 Ackley 3D", _c05)
check("05 Ackley 3D", _val(_c05(), 1e4))

bench("06 Rastrigin 5D", _c06)
check("06 Rastrigin 5D", _val(_c06(), 1e4))

bench("07 Rastrigin 8D", _c07)
check("07 Rastrigin 8D", _val(_c07(), 1e4))

bench("08 Rosenbrock 10D", _c08)
check("08 Rosenbrock 10D", _val(_c08(), 1e4))

bench("09 Schwefel 10D", _c09)
check("09 Schwefel 10D", _val(_c09(), 1e4))
