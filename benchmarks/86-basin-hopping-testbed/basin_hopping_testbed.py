#!/usr/bin/env python3
"""Experiment 86 -- Basin Hopping across the NMinimize test corpus (scipy column).

Same 6 problems as ``basin_hopping_testbed.m``, same labels, same order, solved
with scipy.optimize.basinhopping on the identical bounded box, with a MATCHED
parameterisation: K = 8 independent multi-start runs (kept at the min, mirroring
Mathilda's "SearchPoints" -> 8), and a box-scaled stepsize for the wide-box
deceptive functions (the uniform 0.5 displacement cannot traverse a [-500,500]
box in 100 hops). Each run uses the bounded L-BFGS-B quench.

This is the cleanly-comparable subset of the NMinimize test corpus -- the
continuous, box-bounded problems where both implementations, with the matched
parameterisation, reach the same global. The rest of the corpus is analysed in
README.md: the funnel cases where Mathilda's walk reaches the global and scipy's
stalls even at K=8 (Rastrigin-5D/8D), the deceptive high-D Schwefel both miss,
the constrained / equality / disjunctive cases scipy's box-only basinhopping
cannot express, and the combinatorial / mixed-integer cases.

The check is the objective at the global optimum, rounded per problem so both
implementations, reaching the same basin, agree.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import basinhopping

from harness import bench, check, require

require(["scipy.optimize:basinhopping"])

K = 8


def _solve(fun, bounds, stepsize=0.5):
    # K independent basinhopping runs from seeded random starts in the box, kept
    # at the min -- the multi-start analogue of Mathilda's "SearchPoints" -> K.
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])
    best = None
    for s in range(1, K + 1):
        rng = np.random.default_rng(s)
        x0 = rng.uniform(lo, hi)
        mk = {"method": "L-BFGS-B", "bounds": bounds}
        r = basinhopping(fun, x0, minimizer_kwargs=mk, niter=100, seed=s, stepsize=stepsize)
        if best is None or r.fun < best.fun:
            best = r
    return best


def _val(res, scale):
    return int(np.floor(scale * res.fun + 0.5))


def _chained(v):
    x = v[0]
    return np.sin(2 * x) + np.cos(x)

def _rugged(v):
    x, y = v
    return x * x + y * y + 10 * np.sin(3 * x) ** 2 + 10 * np.sin(3 * y) ** 2

def _schwefel2(v):
    x, y = v
    return 837.9658 - x * np.sin(np.sqrt(abs(x))) - y * np.sin(np.sqrt(abs(y)))

def _schaffer2(v):
    x, y = v
    return 0.5 + (np.sin(x * x - y * y) ** 2 - 0.5) / (1 + 0.001 * (x * x + y * y)) ** 2

def _ackley(v):
    v = np.asarray(v)
    n = v.size
    return (-20 * np.exp(-0.2 * np.sqrt(np.sum(v * v) / n))
            - np.exp(np.sum(np.cos(2 * np.pi * v)) / n) + 20 + np.e)

def _rosenbrock(v):
    v = np.asarray(v)
    return float(np.sum(100 * (v[1:] - v[:-1] ** 2) ** 2 + (1 - v[:-1]) ** 2))


_c01 = lambda: _solve(_chained,    [(-2, 3)])
_c02 = lambda: _solve(_rugged,     [(-5, 5)] * 2)
_c03 = lambda: _solve(_schwefel2,  [(-500, 500)] * 2, stepsize=150)
_c04 = lambda: _solve(_schaffer2,  [(-100, 100)] * 2, stepsize=40)
_c05 = lambda: _solve(_ackley,     [(-32, 32)] * 3)
_c06 = lambda: _solve(_rosenbrock, [(-5, 5)] * 10)


bench("01 Chained trig 1D", _c01)
check("01 Chained trig 1D", _val(_c01(), 1e4))

bench("02 Rugged sine 2D", _c02)
check("02 Rugged sine 2D", _val(_c02(), 1e4))

bench("03 Schwefel 2D", _c03)
check("03 Schwefel 2D", _val(_c03(), 1e3))

bench("04 Schaffer N2 2D", _c04)
check("04 Schaffer N2 2D", _val(_c04(), 1e4))

bench("05 Ackley 3D", _c05)
check("05 Ackley 3D", _val(_c05(), 1e4))

bench("06 Rosenbrock 10D", _c06)
check("06 Rosenbrock 10D", _val(_c06(), 1e4))
