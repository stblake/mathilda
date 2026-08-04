#!/usr/bin/env python3
"""Experiment 13 -- Numerical ODE integration (scipy.integrate.solve_ivp column).

Same five kernels as ``ode_integration.m``, same order.

NON-STIFF AND STIFF ARE DIFFERENT MEASUREMENTS.  On a stiff system an explicit
method is not slow, it is unusable -- the step size collapses to keep stability.
The Van der Pol row therefore uses LSODA (which switches automatically), and its
check is what proves the answer is real rather than merely fast.

Tolerances are set tight (rtol/atol 1e-10) because the checks compare four
decimal places against the CAS columns' default accuracy goals; a loose
tolerance would show up as a check disagreement rather than as a timing.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import numpy as np
from scipy.integrate import solve_ivp

from harness import bench, check, require

require(["scipy.integrate:solve_ivp"])

TOL = dict(rtol=1e-10, atol=1e-12)


def r4(x):
    return int(math.floor(x * 1e4 + 0.5))


def run(f, y0, t1, method="RK45"):
    return solve_ivp(f, (0.0, t1), y0, method=method, **TOL)


bench("NDSolve y'=-y on [0,10]",
      lambda: run(lambda t, y: [-y[0]], [1.0], 10.0))
check("NDSolve y'=-y on [0,10]",
      r4(run(lambda t, y: [-y[0]], [1.0], 10.0).y[0, -1]))

bench("NDSolve harmonic oscillator",
      lambda: run(lambda t, y: [y[1], -y[0]], [1.0, 0.0], 20.0))
check("NDSolve harmonic oscillator",
      r4(run(lambda t, y: [y[1], -y[0]], [1.0, 0.0], 20.0).y[0, -1]))


def vdp(t, y):
    return [y[1], 10.0 * (1.0 - y[0] ** 2) * y[1] - y[0]]


bench("NDSolve Van der Pol mu=10",
      lambda: run(vdp, [2.0, 0.0], 20.0, method="LSODA"), reps=1)
check("NDSolve Van der Pol mu=10",
      r4(run(vdp, [2.0, 0.0], 20.0, method="LSODA").y[0, -1]))

bench("NDSolve y'=y^2 t nonlinear",
      lambda: run(lambda t, y: [y[0] ** 2 * t], [0.5], 1.0))
check("NDSolve y'=y^2 t nonlinear",
      r4(run(lambda t, y: [y[0] ** 2 * t], [0.5], 1.0).y[0, -1]))

bench("NDSolve long horizon, t to 200",
      lambda: run(lambda t, y: [-y[0] / 10.0], [1.0], 200.0))
check("NDSolve long horizon, t to 200",
      r4(run(lambda t, y: [-y[0] / 10.0], [1.0], 200.0).y[0, -1]))
