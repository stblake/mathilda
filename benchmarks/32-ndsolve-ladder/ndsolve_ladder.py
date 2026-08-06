#!/usr/bin/env python3
"""Experiment 32 -- NDSolve difficulty ladder (scipy.integrate column).

Same rows as ``ndsolve_ladder.m``, same order, same labels.

METHOD CHOICE IS PART OF THE MEASUREMENT.  The non-stiff rows use RK45, an
explicit Runge-Kutta of the same family the CAS columns reach for by default.
The Van der Pol rows use LSODA, which switches to BDF when it detects
stiffness -- the same thing a CAS is claiming to do when its mu=1000 row does
not blow up. Handing scipy a stiff solver where the problem is stiff is the
honest baseline: the alternative is timing RK45 failing, which measures
nothing.

TOLERANCES.  rtol=1e-10, atol=1e-12, tight enough that the four-decimal checks
compare solver quality rather than tolerance settings.

WHY EVERY ROW EVALUATES AT THE FINAL TIME.  solve_ivp integrates to t1 or
raises; it cannot quietly return a trajectory that stops short. The CAS columns
can, and the endpoint check is what catches it -- so this column's job on those
rows is to supply the value the CAS should have produced.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import numpy as np
from scipy.integrate import solve_ivp

from harness import bench, check, require

require(["scipy.integrate:solve_ivp"])

TOL = dict(rtol=1e-10, atol=1e-12)


def r(x, p):
    v = float(x) * (10 ** p)
    return int(math.floor(v + 0.5)) if v >= 0 else -int(math.floor(-v + 0.5))


def run(f, y0, t1, method="RK45", t_eval_at=None):
    kw = dict(TOL)
    if t_eval_at is not None:
        kw["t_eval"] = [t_eval_at]
    return solve_ivp(f, (0.0, t1), y0, method=method, **kw)


# ---- horizon, decaying linear -------------------------------------------
def decay(t1, at):
    s = run(lambda t, y: [-y[0]], [1.0], t1, t_eval_at=at)
    return s.y[0, -1]


bench("NDS decay t=10", lambda: run(lambda t, y: [-y[0]], [1.0], 10.0))
check("NDS decay t=10", r(decay(10.0, 10.0), 6))

bench("NDS decay t=1000", lambda: run(lambda t, y: [-y[0]], [1.0], 1000.0))
check("NDS decay t=1000", r(decay(1000.0, 20.0), 6))


# ---- horizon, oscillatory ------------------------------------------------
def harmonic(t1):
    s = run(lambda t, y: [y[1], -y[0]], [1.0, 0.0], t1, t_eval_at=t1)
    return s.y[0, -1]


bench("NDS harmonic t=100",
      lambda: run(lambda t, y: [y[1], -y[0]], [1.0, 0.0], 100.0))
check("NDS harmonic t=100", r(harmonic(100.0), 4))

bench("NDS harmonic t=1000",
      lambda: run(lambda t, y: [y[1], -y[0]], [1.0, 0.0], 1000.0))
check("NDS harmonic t=1000", r(harmonic(1000.0), 4))


# ---- stiffness: Van der Pol ---------------------------------------------
def vdp_f(mu):
    return lambda t, y: [y[1], mu * (1.0 - y[0] ** 2) * y[1] - y[0]]


def vdp(mu):
    s = run(vdp_f(mu), [2.0, 0.0], 20.0, method="LSODA", t_eval_at=20.0)
    return s.y[0, -1]


bench("NDS Van der Pol mu=10",
      lambda: run(vdp_f(10.0), [2.0, 0.0], 20.0, method="LSODA"))
check("NDS Van der Pol mu=10", r(vdp(10.0), 3))

bench("NDS Van der Pol mu=100",
      lambda: run(vdp_f(100.0), [2.0, 0.0], 20.0, method="LSODA"), reps=1)
check("NDS Van der Pol mu=100", r(vdp(100.0), 2))

bench("NDS Van der Pol mu=1000",
      lambda: run(vdp_f(1000.0), [2.0, 0.0], 20.0, method="LSODA"), reps=1)
check("NDS Van der Pol mu=1000", r(vdp(1000.0), 1))


# ---- sensitivity: Lorenz -------------------------------------------------
def lorenz_f(t, y):
    return [10.0 * (y[1] - y[0]),
            y[0] * (28.0 - y[2]) - y[1],
            y[0] * y[1] - (8.0 / 3.0) * y[2]]


def lorenz():
    s = run(lorenz_f, [1.0, 1.0, 1.0], 5.0, t_eval_at=5.0)
    return s.y[0, -1]


bench("NDS Lorenz t=5", lambda: run(lorenz_f, [1.0, 1.0, 1.0], 5.0))
check("NDS Lorenz t=5", r(lorenz(), 3))
