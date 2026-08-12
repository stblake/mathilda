#!/usr/bin/env python3
"""Experiment 12 -- Numerical quadrature (scipy.integrate column).

Same six kernels as ``quadrature.m``, same order.

FOUR INTEGRAND CLASSES, because one adaptive rule cannot be good at all four and
a single smooth row would hide that: smooth, endpoint-singular, oscillatory and
2-D, plus a peaked integrand and an infinite interval.  All have known answers,
so the checks test accuracy as well as agreement.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import numpy as np
from scipy.integrate import dblquad, quad

from harness import bench, check, require

require(["scipy.integrate:quad", "scipy.integrate:dblquad"])


def r(x, p):
    return int(math.floor(x * 10 ** p + 0.5))


bench("NIntegrate Sin[x] on [0,Pi]", lambda: quad(math.sin, 0.0, math.pi))
check("NIntegrate Sin[x] on [0,Pi]", r(quad(math.sin, 0.0, math.pi)[0], 6))

f2 = lambda t: 1.0 / math.sqrt(t) if t > 0 else 0.0
bench("NIntegrate 1/Sqrt[x] on [0,1]", lambda: quad(f2, 0.0, 1.0))
check("NIntegrate 1/Sqrt[x] on [0,1]", r(quad(f2, 0.0, 1.0)[0], 4))

f3 = lambda t: math.sin(40.0 * t)
bench("NIntegrate Sin[40 x] on [0,Pi]", lambda: quad(f3, 0.0, math.pi))
check("NIntegrate Sin[40 x] on [0,Pi]", r(quad(f3, 0.0, math.pi)[0], 4))

f4 = lambda t: math.exp(-1000.0 * (t - 0.5) ** 2)
bench("NIntegrate narrow Gaussian", lambda: quad(f4, 0.0, 1.0))
check("NIntegrate narrow Gaussian", r(quad(f4, 0.0, 1.0)[0], 6))

bench("NIntegrate 2-D smooth",
      lambda: dblquad(lambda y, x: math.sin(x) * math.cos(y), 0.0, math.pi,
                      lambda _: 0.0, lambda _: math.pi / 2))
check("NIntegrate 2-D smooth",
      r(dblquad(lambda y, x: math.sin(x) * math.cos(y), 0.0, math.pi,
                lambda _: 0.0, lambda _: math.pi / 2)[0], 4))

f6 = lambda t: math.exp(-t * t)
bench("NIntegrate Exp[-x^2] on [0,Infinity]", lambda: quad(f6, 0.0, np.inf))
check("NIntegrate Exp[-x^2] on [0,Infinity]", r(quad(f6, 0.0, np.inf)[0], 4))
