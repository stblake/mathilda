#!/usr/bin/env python3
"""Experiment 80 -- SHGO across the NMinimize test corpus (scipy column).

Same 11 problems as ``shgo_testbed.m``, same labels, same order, solved with
scipy.optimize.shgo using the SAME sampling_method and sample count (n) and an
equivalent bounded box + constraints, so the race is SHGO-vs-shgo. This is the
cleanly-comparable set (the multimodal stress functions and the
inequality-constrained problems). The rest of the corpus is analysed in
README.md -- the constrained problems scipy.optimize.shgo mishandles, the
high-dimensional / sharp-needle functions outside SHGO's envelope, and the
combinatorial-integer problems (scipy.optimize.shgo is continuous-only).

The check is the objective at the global optimum, rounded per problem so that
both implementations, reaching the same basin, agree.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import shgo

from harness import bench, check, require

require(["scipy.optimize:shgo"])


def _solve(fun, bounds, n, method, constraints=()):
    return shgo(fun, bounds, n=n, iters=1, sampling_method=method,
                constraints=constraints)


def _val(res, scale):
    f = res.fun
    if f is None or not np.isfinite(f):
        return "NOFEASIBLE"
    return int(np.floor(scale * f + 0.5))


def _eggholder(v):
    x, y = v
    return (-(y + 47) * np.sin(np.sqrt(abs(y + x / 2 + 47)))
            - x * np.sin(np.sqrt(abs(x - (y + 47)))))

def _schwefel2(v):
    x, y = v
    return 837.9658 - x * np.sin(np.sqrt(abs(x))) - y * np.sin(np.sqrt(abs(y)))

def _schaffer2(v):
    x, y = v
    return 0.5 + (np.sin(x * x - y * y) ** 2 - 0.5) / (1 + 0.001 * (x * x + y * y)) ** 2

def _mishra(v):
    x, y = v
    return (np.sin(y) * np.exp((1 - np.cos(x)) ** 2)
            + np.cos(x) * np.exp((1 - np.sin(y)) ** 2) + (x - y) ** 2)

def _rugged(v):
    x, y = v
    return x * x + y * y + 10 * np.sin(3 * x) ** 2 + 10 * np.sin(3 * y) ** 2

def _ackley3(v):
    v = np.asarray(v)
    return (-20 * np.exp(-0.2 * np.sqrt(np.sum(v * v) / 3))
            - np.exp(np.sum(np.cos(2 * np.pi * v)) / 3) + 20 + np.e)

def _rastrigin5(v):
    v = np.asarray(v)
    return 50 + float(np.sum(v * v - 10 * np.cos(2 * np.pi * v)))

def _griewank5(v):
    v = np.asarray(v)
    idx = np.arange(1, 6)
    return float(np.sum(v * v) / 4000 - np.prod(np.cos(v / np.sqrt(idx))) + 1)


_c01 = lambda: _solve(_eggholder,  [(-512, 512)] * 2,   500, "simplicial")
_c02 = lambda: _solve(_schwefel2,  [(-500, 500)] * 2,   300, "simplicial")
_c03 = lambda: _solve(_schaffer2,  [(-100, 100)] * 2,   400, "simplicial")
_c04 = lambda: _solve(_mishra,     [(-10, 10)] * 2,     300, "simplicial",
                      ({"type": "ineq", "fun": lambda v: 25 - (v[0] + 5) ** 2 - (v[1] + 5) ** 2},))
_c05 = lambda: _solve(_rugged,     [(-5, 5)] * 2,       200, "simplicial")
_c06 = lambda: _solve(_ackley3,    [(-32, 32)] * 3,     200, "simplicial")
_c07 = lambda: _solve(_rastrigin5, [(-5.12, 5.12)] * 5, 300, "simplicial")
_c08 = lambda: _solve(_griewank5,  [(-15, 15)] * 5,     300, "simplicial")
_c09 = lambda: _solve(lambda v: v[0] + v[1], [(-10, 10)] * 2, 100, "simplicial",
                      ({"type": "ineq", "fun": lambda v: 9 - v[0] ** 2 - v[1] ** 2},))
_c10 = lambda: _solve(lambda v: np.sin(2 * v[0]) + np.cos(v[0]), [(-2, 3)], 100, "simplicial")
_c11 = lambda: _solve(lambda v: v[0] ** 2 + v[1] ** 2, [(-10, 10)] * 2, 100, "simplicial",
                      ({"type": "ineq", "fun": lambda v: v[0] + v[1] - 1},))

bench("01 Eggholder 2D", _c01)
check("01 Eggholder 2D", _val(_c01(), 1e0))

bench("02 Schwefel 2D", _c02)
check("02 Schwefel 2D", _val(_c02(), 1e3))

bench("03 Schaffer N2 2D", _c03)
check("03 Schaffer N2 2D", _val(_c03(), 1e4))

bench("04 Mishra Bird disk", _c04)
check("04 Mishra Bird disk", _val(_c04(), 1e3))

bench("05 Rugged sine 2D", _c05)
check("05 Rugged sine 2D", _val(_c05(), 1e4))

bench("06 Ackley 3D", _c06)
check("06 Ackley 3D", _val(_c06(), 1e4))

bench("07 Rastrigin 5D", _c07)
check("07 Rastrigin 5D", _val(_c07(), 1e4))

bench("08 Griewank 5D", _c08)
check("08 Griewank 5D", _val(_c08(), 1e4))

bench("09 Disk linear", _c09)
check("09 Disk linear", _val(_c09(), 1e4))

bench("10 Chained trig 1D", _c10)
check("10 Chained trig 1D", _val(_c10(), 1e4))

bench("11 Quadratic halfplane", _c11)
check("11 Quadratic halfplane", _val(_c11(), 1e4))
