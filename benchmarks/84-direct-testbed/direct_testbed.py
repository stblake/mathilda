#!/usr/bin/env python3
"""Experiment 84 -- DIRECT across the NMinimize test corpus (scipy column).

Same 18 problems as ``direct_testbed.m``, same labels, same order, solved with
scipy.optimize.direct on the identical bounded box. DIRECT is deterministic --
no seed. This is the RACEABLE subset of the NMinimize test corpus
(tests/test_nminimize.c): the continuous, box-bounded problems where the two
DIRECT implementations do the same work. Both columns use scipy's default
parameters except where a case is annotated ``LF/N`` (locally_biased=False,
maxfun=N) -- the same appropriate parameterisation on both sides.

The race is raw-vs-raw (scipy's direct does no local polish, so the .m column
sets "PostProcess" -> False). Every case reaches the SAME point in both systems:
the global for most, a shared non-global basin for Rosenbrock-10D (~8.7916) and
Gaussian-well-10D (both raw DIRECTs stall at 1.0 -- the well at 1.2345 is too
narrow for center-based subdivision).

The rest of the corpus is analysed in README.md and is NOT a DIRECT-vs-direct
race: scipy.optimize.direct is BOX-ONLY, so the constrained / equality /
disjunctive / integer problems have no scipy column (Mathilda's DIRECT solves
the low-dimensional ones via its penalty + augmented-Lagrangian polish, which
scipy's direct cannot do at all); the high-dimensional equality-constrained
A-series and the combinatorial B-series are beyond DIRECT in either system.

The check is the objective at the reported point, rounded per problem so both
implementations agree at the basin.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import direct
import scipy.special as _sp

from harness import bench, check, require

require(["scipy.optimize:direct"])


def _solve(fun, bounds, **kw):
    return direct(fun, bounds, **kw)


def _val(res, scale):
    return int(np.floor(scale * res.fun + 0.5))


# ---- objectives ----------------------------------------------------------

def _chained(v):    x = v[0]; return np.sin(2 * x) + np.cos(x)
def _gamma(v):      return float(_sp.gamma(v[0]))
def _quartic(v):    x = v[0]; return x ** 4 - 3 * x ** 2 - x

def _schwefel2(v):
    x, y = v
    return 837.9658 - x * np.sin(np.sqrt(abs(x))) - y * np.sin(np.sqrt(abs(y)))

def _schaffer2(v):
    x, y = v
    return 0.5 + (np.sin(x * x - y * y) ** 2 - 0.5) / (1 + 0.001 * (x * x + y * y)) ** 2

def _rugged(v):
    x, y = v
    return x * x + y * y + 10 * np.sin(3 * x) ** 2 + 10 * np.sin(3 * y) ** 2

def _eggholder(v):
    x, y = v
    return (-(y + 47) * np.sin(np.sqrt(abs(y + x / 2 + 47)))
            - x * np.sin(np.sqrt(abs(x - (y + 47)))))

def _griewank(v):
    v = np.asarray(v)
    return (np.sum(v ** 2 / 4000.0)
            - np.prod(np.cos(v / np.sqrt(np.arange(1, v.size + 1)))) + 1.0)

def _scaledquad(v):
    return float(sum(10.0 ** (6 * i / 9.0) * v[i] ** 2 for i in range(len(v))))

def _quadbowl(v):
    v = np.asarray(v); return float(np.sum(v * v))

def _ackley(v):
    v = np.asarray(v); n = v.size
    return (-20.0 * np.exp(-0.2 * np.sqrt(np.sum(v * v) / n))
            - np.exp(np.sum(np.cos(2 * np.pi * v)) / n) + 20.0 + np.e)

def _katsuura(v):
    v = np.asarray(v); d = v.size
    prod = 1.0
    for i in range(d):
        s = 1.0 + (i + 1) * sum(abs(2 ** k * v[i] - round(2 ** k * v[i])) / 2 ** k
                                for k in range(1, 26))
        prod *= s ** (10.0 / d ** 1.2)
    return (10.0 / d ** 2) * prod - (10.0 / d ** 2)

def _rastrigin(v):
    v = np.asarray(v)
    return 10.0 * v.size + np.sum(v * v - 10.0 * np.cos(2 * np.pi * v))

def _rosenbrock(v):
    v = np.asarray(v)
    return float(np.sum(100 * (v[1:] - v[:-1] ** 2) ** 2 + (1 - v[:-1]) ** 2))

def _modackley(v):
    v = np.asarray(v); d = v.size
    return (-np.exp(-0.2 * np.sqrt(np.sum(v * v) / d)) * np.prod(np.cos(20 * v))
            + 0.05 * np.sum(v * v))

def _gwell(v):
    v = np.asarray(v)
    return 1.0 - np.exp(-2.0 * np.sum((v - 1.2345) ** 2)) + 1e-5 * np.sum(v ** 2)


# ---- cases ---------------------------------------------------------------

_c01 = lambda: _solve(_chained,    [(-2, 3)])
_c02 = lambda: _solve(_quartic,    [(-5, 5)])
_c03 = lambda: _solve(_gamma,      [(1, 2)])
_c04 = lambda: _solve(_schwefel2,  [(-500, 500)] * 2)
_c05 = lambda: _solve(_schaffer2,  [(-100, 100)] * 2)
_c06 = lambda: _solve(_rugged,     [(-5, 5)] * 2)
_c07 = lambda: _solve(_eggholder,  [(-512, 512)] * 2, locally_biased=False, maxfun=20000)
_c08 = lambda: _solve(_griewank,   [(-15, 15)] * 5)
_c09 = lambda: _solve(_griewank,   [(-600, 600)] * 10)
_c10 = lambda: _solve(_scaledquad, [(-10, 10)] * 10)
_c11 = lambda: _solve(_quadbowl,   [(-5, 5)] * 3)
_c12 = lambda: _solve(_ackley,     [(-32, 32)] * 3)
_c13 = lambda: _solve(_katsuura,   [(-100, 100)] * 8)
_c14 = lambda: _solve(_rastrigin,  [(-5.12, 5.12)] * 5)
_c15 = lambda: _solve(_rastrigin,  [(-5.12, 5.12)] * 8)
_c16 = lambda: _solve(_rosenbrock, [(-5, 5)] * 10)
_c17 = lambda: _solve(_modackley,  [(-5, 5)] * 10)
_c18 = lambda: _solve(_gwell,      [(-5, 5)] * 10, locally_biased=False, maxfun=40000)

bench("01 Chained trig 1D", _c01);         check("01 Chained trig 1D", _val(_c01(), 1e4))
bench("02 Quartic 1D", _c02);              check("02 Quartic 1D", _val(_c02(), 1e4))
bench("03 Gamma 1D", _c03);                check("03 Gamma 1D", _val(_c03(), 1e5))
bench("04 Schwefel 2D", _c04);             check("04 Schwefel 2D", _val(_c04(), 1e3))
bench("05 Schaffer N2 2D", _c05);          check("05 Schaffer N2 2D", _val(_c05(), 1e4))
bench("06 Rugged sine 2D", _c06);          check("06 Rugged sine 2D", _val(_c06(), 1e4))
bench("07 Eggholder 2D LF/20000", _c07);   check("07 Eggholder 2D LF/20000", _val(_c07(), 1e3))
bench("08 Griewank 5D", _c08);             check("08 Griewank 5D", _val(_c08(), 1e4))
bench("09 Griewank 10D", _c09);            check("09 Griewank 10D", _val(_c09(), 1e4))
bench("10 Scaled-quadratic 10D", _c10);    check("10 Scaled-quadratic 10D", _val(_c10(), 1e4))
bench("11 Quadratic bowl 3D", _c11);       check("11 Quadratic bowl 3D", _val(_c11(), 1e4))
bench("12 Ackley 3D", _c12);               check("12 Ackley 3D", _val(_c12(), 1e4))
bench("13 Katsuura 8D", _c13);             check("13 Katsuura 8D", _val(_c13(), 1e4))
bench("14 Rastrigin 5D", _c14);            check("14 Rastrigin 5D", _val(_c14(), 1e4))
bench("15 Rastrigin 8D", _c15);            check("15 Rastrigin 8D", _val(_c15(), 1e4))
bench("16 Rosenbrock 10D", _c16);          check("16 Rosenbrock 10D", _val(_c16(), 1e3))
bench("17 Modified Ackley 10D", _c17);     check("17 Modified Ackley 10D", _val(_c17(), 1e4))
bench("18 Gaussian well 10D LF/40000", _c18); check("18 Gaussian well 10D LF/40000", _val(_c18(), 1e3))
