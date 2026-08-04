#!/usr/bin/env python3
"""Experiment 16 -- Interpolation (scipy.interpolate column).

Same five kernels as ``interpolation.m``, same order.

BUILD AND EVALUATE ARE SEPARATE ROWS on purpose: a system can be fast to
construct and slow to evaluate, and the two live in different files.

Mathematica's Interpolation defaults to order 3 (cubic), so CubicSpline is the
matching constructor.  The two use different end conditions (not-a-knot vs
Wolfram's), which is why the checks evaluate well inside the interval rather
than near an endpoint.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import numpy as np
from scipy.interpolate import CubicSpline, RegularGridInterpolator, interp1d

from harness import bench, check, require, seed

require(["scipy.interpolate:CubicSpline", "scipy.interpolate:interp1d",
         "scipy.interpolate:RegularGridInterpolator"])

xi = np.arange(1, 2001) / 10.0
yi = np.sin(xi)


def r6(x):
    return int(math.floor(float(x) * 1e6 + 0.5))


bench("Interpolation build, 2000 knots", lambda: CubicSpline(xi, yi))
check("Interpolation build, 2000 knots", len(xi))

cs = CubicSpline(xi, yi)
ev = np.arange(1000, 21000) / 1000.0
bench("Interpolation evaluate, 20000 points", lambda: cs(ev))
check("Interpolation evaluate, 20000 points", r6(cs(2.5)))

bench("Interpolation order 1 (linear) build",
      lambda: interp1d(xi, yi, kind="linear"))
check("Interpolation order 1 (linear) build",
      r6(interp1d(xi, yi, kind="linear")(2.5)))

gi = np.arange(1, 61) / 10.0
grid = np.outer(np.sin(gi), np.cos(gi))
ax = np.arange(1.0, 61.0)
bench("ListInterpolation 2-D 60x60",
      lambda: RegularGridInterpolator((ax, ax), grid, method="cubic"))
_rg = RegularGridInterpolator((ax, ax), grid, method="cubic")
check("ListInterpolation 2-D 60x60", r6(_rg([[5.0, 5.0]])[0]))

seed()
xs = np.random.random(100000) * 100.0 + 1.0
bench("Interpolation over 10^5 array", lambda: cs(xs))
check("Interpolation over 10^5 array", r6(cs(7.5)))
