#!/usr/bin/env python3
"""Experiment 16 -- Interpolation (scipy.interpolate column).

Same kernels as ``interpolation.m``, same labels: the ``Interpolation`` rows,
then the value-only ``ListInterpolation`` companion (1-D build, 1-D vectorised
evaluate, 2-D build).

BUILD AND EVALUATE ARE SEPARATE ROWS on purpose: a system can be fast to
construct and slow to evaluate, and the two live in different files.

Mathematica's Interpolation defaults to order 3 (cubic), so CubicSpline is the
matching constructor; ListInterpolation over a value array is the same cubic on
a synthesised grid, so CubicSpline over an ``arange`` grid is its analog, and
RegularGridInterpolator(method="cubic") for the tensor case.  The two use
different end conditions (not-a-knot vs Wolfram's), so EVERY check evaluates at a
GRID NODE -- a data point, where every scheme returns the stored value -- rather
than between nodes, where the schemes legitimately disagree.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline, interp1d

from harness import bench, check, require, seed

# RectBivariateSpline, not RegularGridInterpolator(method="cubic"), is the tensor
# analog: both Mathilda's 2-D ListInterpolation and RectBivariateSpline are
# node-EXACT interpolating tensor cubics (RectBivariateSpline to 3e-16 at every
# node), so a grid-node check agrees.  RegularGridInterpolator's cubic method is
# a smoothing tensor spline that misses the data value at a node by up to ~4e-5,
# which silently CHECK-FAIL'd the 2-D row and discarded its timing.
require(["scipy.interpolate:CubicSpline", "scipy.interpolate:interp1d",
         "scipy.interpolate:RectBivariateSpline"])

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
      lambda: RectBivariateSpline(ax, ax, grid))
_rg = RectBivariateSpline(ax, ax, grid)
check("ListInterpolation 2-D 60x60", r6(_rg.ev(5.0, 5.0)))   # node -> sin(0.5)cos(0.5)

seed()
xs = np.random.random(100000) * 100.0 + 1.0
bench("Interpolation over 10^5 array", lambda: cs(xs))
check("Interpolation over 10^5 array", r6(cs(7.5)))

# -------------------------------------------------------------------------
# ListInterpolation -- values on a regular grid.  Mathilda synthesises the
# integer grid 1, 2, ... and delegates to the cubic Interpolation engine, so
# CubicSpline over np.arange(1, n+1) is the 1-D analog and
# RegularGridInterpolator(method="cubic") the tensor analog.
# -------------------------------------------------------------------------

# 1-D build from 10^5 values on the integer grid.
lxi = np.arange(1, 100001)
lvals = np.sin(lxi / 10.0)
bench("ListInterpolation 1-D build, 10^5 values", lambda: CubicSpline(lxi, lvals))
lcs = CubicSpline(lxi, lvals)
check("ListInterpolation 1-D build, 10^5 values", r6(lcs(25)))   # node 25 -> sin(2.5)

# 1-D vectorised evaluate over a 10^5 in-domain query array.
seed()
lxs = np.random.random(100000) * 99000.0 + 1.0
bench("ListInterpolation 1-D evaluate, 10^5 array", lambda: lcs(lxs))
check("ListInterpolation 1-D evaluate, 10^5 array", r6(lcs(7500)))  # node -> sin(750)

# 2-D build: a 200x200 separable grid.
lgi = np.arange(1, 201) / 10.0
lgrid = np.outer(np.sin(lgi), np.cos(lgi))
lax = np.arange(1.0, 201.0)
bench("ListInterpolation 2-D build, 200x200",
      lambda: RectBivariateSpline(lax, lax, lgrid))
_lrg = RectBivariateSpline(lax, lax, lgrid)
check("ListInterpolation 2-D build, 200x200", r6(_lrg.ev(50.0, 80.0)))  # node -> sin(5)cos(8)
