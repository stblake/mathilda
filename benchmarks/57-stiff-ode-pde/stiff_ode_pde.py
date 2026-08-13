#!/usr/bin/env python3
"""Experiment 57 -- Stiff ODEs and PDEs (scipy.integrate column).

Same kernels as ``stiff_ode_pde.m``, same order and labels.

Baseline is compiled (scipy.integrate.solve_ivp / LSODA / BDF over
ODEPACK/Fortran).  The PDE rows use a hand-written method-of-lines discretisation
(second-difference stencil + solve_ivp), matching NDSolve's MethodOfLines.

Checks are integrator-robust (steady-state / analytic values); the heat and wave
rows check the exact separable solutions exp(-pi^2 t) sin(pi x) and
cos(pi t) sin(pi x) at (0.1, 0.5), which any accurate solver reproduces.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math
import numpy as np
from scipy.integrate import solve_ivp

from harness import bench, check, require

require(["scipy.integrate:solve_ivp"])


def r(x, p):
    return int(math.floor(float(x) * 10 ** p + 0.5))


# ---- stiff / non-stiff ODEs ---------------------------------------------

def stiff_scalar():
    return solve_ivp(lambda t, y: -1000 * y + 3000 - 2000 * math.exp(-t),
                     (0, 1), [0.0], method="BDF", rtol=1e-8, atol=1e-10, t_eval=[1.0])

bench("NDSolve stiff scalar relaxation", stiff_scalar)
check("NDSolve stiff scalar relaxation", r(stiff_scalar().y[0, -1], 4))


def robertson():
    def f(t, y):
        a, b, c = y
        return [-0.04 * a + 1e4 * b * c,
                0.04 * a - 1e4 * b * c - 3e7 * b * b,
                3e7 * b * b]
    return solve_ivp(f, (0, 1), [1.0, 0.0, 0.0], method="Radau",
                     rtol=1e-8, atol=1e-10, t_eval=[1.0])

bench("NDSolve Robertson stiff kinetics", robertson)
check("NDSolve Robertson stiff kinetics", r(robertson().y[2, -1], 4))


def harmonic():
    return solve_ivp(lambda t, y: [y[1], -y[0]], (0, 1), [1.0, 0.0],
                     method="RK45", rtol=1e-10, atol=1e-12, t_eval=[1.0])

bench("NDSolve linear system harmonic", harmonic)
check("NDSolve linear system harmonic", r(harmonic().y[0, -1], 6))


def vdp():
    return solve_ivp(lambda t, y: [y[1], 10 * (1 - y[0] ** 2) * y[1] - y[0]],
                     (0, 1), [2.0, 0.0], method="LSODA", rtol=1e-9, atol=1e-11,
                     t_eval=[1.0])

bench("NDSolve Van der Pol mu=10", vdp)
check("NDSolve Van der Pol mu=10", r(vdp().y[0, -1], 3))


# ---- PDEs via method of lines -------------------------------------------

N = 100
xs = np.linspace(0.0, 1.0, N + 1)          # 0, 0.01, ..., 1; x=0.5 at index 50
dx = 1.0 / N
xi = xs[1:-1]                               # interior nodes


def d2(u):                                  # second difference, Dirichlet BC (u=0 at ends)
    lap = np.empty_like(u)
    lap[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / dx ** 2
    lap[0] = (-2 * u[0] + u[1]) / dx ** 2
    lap[-1] = (u[-2] - 2 * u[-1]) / dx ** 2
    return lap


def heat():
    u0 = np.sin(np.pi * xi)
    return solve_ivp(lambda t, u: d2(u), (0, 0.1), u0, method="BDF",
                     rtol=1e-8, atol=1e-10, t_eval=[0.1])

bench("NDSolve heat PDE (method of lines)", heat)
check("NDSolve heat PDE (method of lines)", r(heat().y[N // 2 - 1, -1], 3))


def wave():
    u0 = np.sin(np.pi * xi)
    y0 = np.concatenate([u0, np.zeros_like(u0)])

    def f(t, y):
        u, v = y[:len(xi)], y[len(xi):]
        return np.concatenate([v, d2(u)])
    return solve_ivp(f, (0, 0.1), y0, method="RK45", rtol=1e-9, atol=1e-11,
                     t_eval=[0.1])

bench("NDSolve wave PDE (method of lines)", wave)
check("NDSolve wave PDE (method of lines)", r(wave().y[N // 2 - 1, -1], 4))
