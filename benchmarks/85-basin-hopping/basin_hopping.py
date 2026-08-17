#!/usr/bin/env python3
"""Experiment 85 -- Basin Hopping global optimization (scipy column).

Same 8 problems as ``basin_hopping.m``, same labels, same order, solved with
scipy.optimize.basinhopping using scipy's DEFAULT parameters (T 1, stepsize 0.5,
interval 50, target_accept_rate 0.5, stepwise_factor 0.9, niter 100), a single run
from a seeded random start in the box, and L-BFGS-B as the bounded local minimizer
(the "quench"). The race is BasinHopping-vs-basinhopping on the same algorithm.
These are standard bounded multimodal benchmarks -- Himmelblau, Booth, Beale,
Rosenbrock, six-hump camel, Ackley (2-D and 5-D), sphere -- whose global both
systems' single run reaches and polishes to agreement. The check is the objective
at the global optimum rounded to 10^6, matching the .m ``Round[10^6 First[...]]``.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import basinhopping

from harness import bench, check, require

require(["scipy.optimize:basinhopping"])

SEED = 1


def _solve(fun, bounds):
    # A single basinhopping run from a seeded random start in the box, with the
    # bounded L-BFGS-B quench -- scipy's default niter/T/stepsize on both sides.
    rng = np.random.default_rng(SEED)
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])
    x0 = rng.uniform(lo, hi)
    mk = {"method": "L-BFGS-B", "bounds": bounds}
    return basinhopping(fun, x0, minimizer_kwargs=mk, niter=100, seed=SEED)


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

def _sphere(v):
    x, y = v
    return x * x + y * y


_c01 = lambda: _solve(_himmelblau,  [(-5, 5), (-5, 5)])
_c02 = lambda: _solve(_booth,       [(-10, 10), (-10, 10)])
_c03 = lambda: _solve(_beale,       [(-4.5, 4.5), (-4.5, 4.5)])
_c04 = lambda: _solve(_rosenbrock,  [(-5, 5), (-5, 5)])
_c05 = lambda: _solve(_camel,       [(-3, 3), (-2, 2)])
_c06 = lambda: _solve(_ackley,      [(-5, 5), (-5, 5)])
_c07 = lambda: _solve(_sphere,      [(-5, 5), (-5, 5)])
_c08 = lambda: _solve(_ackley,      [(-5, 5)] * 5)


def _chk(res):
    return int(np.floor(1e6 * res.fun + 0.5))


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

bench("07 Sphere 2D", _c07)
check("07 Sphere 2D", _chk(_c07()))

bench("08 Ackley 5D", _c08)
check("08 Ackley 5D", _chk(_c08()))
