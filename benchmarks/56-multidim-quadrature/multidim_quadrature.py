#!/usr/bin/env python3
"""Experiment 56 -- Multidimensional & hard quadrature (scipy.integrate column).

Same kernels as ``multidim_quadrature.m``, same order and labels.

Baseline is compiled (scipy.integrate), but dblquad/tplquad call the PYTHON
integrand once per node, so the cubature rows are partly compiled-vs-interpreted
integrand -- the 1-D rows (oscillatory, singular, semi-infinite) are the cleaner
compiled-vs-compiled comparison.  Every integrand has a known closed value.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math
import numpy as np
from scipy.integrate import quad, dblquad, tplquad

from harness import bench, check, require

require(["scipy.integrate:quad", "scipy.integrate:dblquad", "scipy.integrate:tplquad"])


def r6(x):
    return int(math.floor(float(x) * 1e6 + 0.5))


bench("NIntegrate 2D x y on unit square",
      lambda: dblquad(lambda y, x: x * y, 0, 1, 0, 1))
check("NIntegrate 2D x y on unit square",
      r6(dblquad(lambda y, x: x * y, 0, 1, 0, 1)[0]))

bench("NIntegrate 3D x y z on unit cube",
      lambda: tplquad(lambda z, y, x: x * y * z, 0, 1, 0, 1, 0, 1))
check("NIntegrate 3D x y z on unit cube",
      r6(tplquad(lambda z, y, x: x * y * z, 0, 1, 0, 1, 0, 1)[0]))

bench("NIntegrate 2D Gaussian on [-3,3]^2",
      lambda: dblquad(lambda y, x: math.exp(-(x * x + y * y)), -3, 3, -3, 3))
check("NIntegrate 2D Gaussian on [-3,3]^2",
      r6(dblquad(lambda y, x: math.exp(-(x * x + y * y)), -3, 3, -3, 3)[0]))

bench("NIntegrate triangle (dependent bounds)",
      lambda: dblquad(lambda y, x: 1.0, 0, 1, 0, lambda x: x))
check("NIntegrate triangle (dependent bounds)",
      r6(dblquad(lambda y, x: 1.0, 0, 1, 0, lambda x: x)[0]))

bench("NIntegrate unit disk (dependent bounds)",
      lambda: dblquad(lambda y, x: 1.0, -1, 1,
                      lambda x: -math.sqrt(1 - x * x), lambda x: math.sqrt(1 - x * x)))
check("NIntegrate unit disk (dependent bounds)",
      r6(dblquad(lambda y, x: 1.0, -1, 1,
                 lambda x: -math.sqrt(1 - x * x), lambda x: math.sqrt(1 - x * x))[0]))

bench("NIntegrate oscillatory Sin[50x]/x on [1,10]",
      lambda: quad(lambda x: math.sin(50 * x) / x, 1, 10))
check("NIntegrate oscillatory Sin[50x]/x on [1,10]",
      r6(quad(lambda x: math.sin(50 * x) / x, 1, 10)[0]))

bench("NIntegrate singular 1/Sqrt[x] on [0,1]",
      lambda: quad(lambda x: 1.0 / math.sqrt(x), 0, 1))
check("NIntegrate singular 1/Sqrt[x] on [0,1]",
      r6(quad(lambda x: 1.0 / math.sqrt(x), 0, 1)[0]))

bench("NIntegrate semi-infinite Exp[-x^2] on [0,Inf]",
      lambda: quad(lambda x: math.exp(-x * x), 0, np.inf))
check("NIntegrate semi-infinite Exp[-x^2] on [0,Inf]",
      r6(quad(lambda x: math.exp(-x * x), 0, np.inf)[0]))
