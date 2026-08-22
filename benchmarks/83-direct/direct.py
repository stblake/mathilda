#!/usr/bin/env python3
"""Experiment 83 -- DIRECT (DIviding RECTangles) global optimization (scipy col).

Same 9 problems as ``direct.m``, same labels, same order, solved with
scipy.optimize.direct using scipy's DEFAULT parameters (locally_biased=True,
eps 1e-4, maxiter 1000, maxfun 1000*n). DIRECT is deterministic -- there is no
seed. The race is DIRECT-vs-direct on the same algorithm, raw on both sides
(scipy's direct does no local polish, so the .m column disables its polish with
"PostProcess" -> False for an apples-to-apples comparison). These are standard
bounded multimodal benchmarks -- Himmelblau, Booth, Beale, Rosenbrock, six-hump
camel, Ackley (2-D and 5-D), Styblinski-Tang, Rastrigin -- whose global basin
both systems reach. The check is the objective at the reported point rounded to
10^3 (coarser than the polished-engine benchmarks: two raw DIRECT runs agree at
the basin, not to machine precision), matching the .m ``Round[10^3 First[...]]``.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import direct

from harness import bench, check, require

require(["scipy.optimize:direct"])


def _solve(fun, bounds):
    return direct(fun, bounds)


def _himmelblau(v):
    x, y = v
    return (x * x + y - 11.0) ** 2 + (x + y * y - 7.0) ** 2

def _booth(v):
    x, y = v
    return (x + 2 * y - 7.0) ** 2 + (2 * x + y - 5.0) ** 2

def _beale(v):
    x, y = v
    return ((1.5 - x + x * y) ** 2 + (2.25 - x + x * y ** 2) ** 2
            + (2.625 - x + x * y ** 3) ** 2)

def _rosenbrock(v):
    x, y = v
    return 100.0 * (y - x * x) ** 2 + (1.0 - x) ** 2

def _camel(v):
    x, y = v
    return (4 - 2.1 * x * x + x ** 4 / 3) * x * x + x * y + (-4 + 4 * y * y) * y * y

def _ackley(v):
    v = np.asarray(v)
    n = v.size
    return (-20.0 * np.exp(-0.2 * np.sqrt(np.sum(v * v) / n))
            - np.exp(np.sum(np.cos(2 * np.pi * v)) / n) + 20.0 + np.e)

def _styblinski(v):
    v = np.asarray(v)
    return float(np.sum(v ** 4 - 16 * v * v + 5 * v) / 2.0)

def _rastrigin(v):
    v = np.asarray(v)
    return 10.0 * v.size + np.sum(v * v - 10.0 * np.cos(2 * np.pi * v))


_c01 = lambda: _solve(_himmelblau,  [(-5, 5), (-5, 5)])
_c02 = lambda: _solve(_booth,       [(-10, 10), (-10, 10)])
_c03 = lambda: _solve(_beale,       [(-4.5, 4.5), (-4.5, 4.5)])
_c04 = lambda: _solve(_rosenbrock,  [(-5, 5), (-5, 5)])
_c05 = lambda: _solve(_camel,       [(-3, 3), (-2, 2)])
_c06 = lambda: _solve(_ackley,      [(-5, 5), (-5, 5)])
_c07 = lambda: _solve(_ackley,      [(-5, 5)] * 5)
_c08 = lambda: _solve(_styblinski,  [(-5, 5), (-5, 5)])
_c09 = lambda: _solve(_rastrigin,   [(-5.12, 5.12), (-5.12, 5.12)])


def _chk(res):
    return int(np.floor(1e3 * res.fun + 0.5))


bench("01 Himmelblau 2D", _c01)
check("01 Himmelblau 2D", _chk(_c01()))

bench("02 Booth 2D", _c02)
check("02 Booth 2D", _chk(_c02()))

bench("03 Beale 2D", _c03)
check("03 Beale 2D", _chk(_c03()))

bench("04 Rosenbrock 2D", _c04)
check("04 Rosenbrock 2D", _chk(_c04()))

bench("05 Six-hump camel 2D", _c05)
check("05 Six-hump camel 2D", _chk(_c05()))

bench("06 Ackley 2D", _c06)
check("06 Ackley 2D", _chk(_c06()))

bench("07 Ackley 5D", _c07)
check("07 Ackley 5D", _chk(_c07()))

bench("08 Styblinski-Tang 2D", _c08)
check("08 Styblinski-Tang 2D", _chk(_c08()))

bench("09 Rastrigin 2D", _c09)
check("09 Rastrigin 2D", _chk(_c09()))
