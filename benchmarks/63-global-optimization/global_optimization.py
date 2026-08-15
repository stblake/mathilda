#!/usr/bin/env python3
"""Experiment 63 -- Global / constrained optimization (scipy.optimize column).

Same cases as ``global_optimization.m``, same labels, same order. The check is
the objective value at the optimum, rounded -- invariant under which feasible
optimum the solver found. The matching scipy engine is chosen by the method
CHARACTER of the Mathilda case: a global search (DifferentialEvolution /
SimulatedAnnealing) races scipy's differential_evolution / dual_annealing; a
local/polish shape races scipy.optimize.minimize.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import differential_evolution, NonlinearConstraint

from harness import bench, check, require

require(["scipy.optimize:differential_evolution", "scipy.optimize:dual_annealing",
         "scipy.optimize:minimize"])


# ---- A1: refinery pooling (bilinear equality) -- scipy differential_evolution ----
# vars: px, py, xA, xB, yA, yB, cA, cB  (indices 0..7). Global optimum -3900.
def _a1_cost(v):
    px, py, xA, xB, yA, yB, cA, cB = v
    return 6 * cA + 16 * cB + 10 * (xA + xB) - 9 * (xA + yA) - 15 * (xB + yB)


def _a1_eqs(v):
    px, py, xA, xB, yA, yB, cA, cB = v
    return [cA + cB - (px + py), 3 * cA + cB - px * xA - py * yA]


def _a1_ineqs(v):  # <= 0
    px, py, xA, xB, yA, yB, cA, cB = v
    return [xA + yA - 100, xB + yB - 200]


_a1_bnds = [(0, 2.5), (0, 1.5), (0, 100), (0, 200), (0, 100), (0, 200), (0, 300), (0, 300)]
_a1_nlc = [NonlinearConstraint(_a1_eqs, 0, 0), NonlinearConstraint(_a1_ineqs, -np.inf, 0)]


def _a1_solve():
    return differential_evolution(_a1_cost, _a1_bnds, constraints=_a1_nlc,
                                  seed=1, tol=1e-7, polish=True, maxiter=200)


bench("A1 refinery pooling (DE)", _a1_solve)
check("A1 refinery pooling (DE)", int(np.floor(_a1_solve().fun + 0.5)))


# ---- A2: isolated Gaussian well in a flat 10-D landscape -- differential_evolution ----
# Global 1.524e-4 at x_i=1.2345; elsewhere the ~1.0 plateau. popsize*d = 40*10 = 400
# matches Mathilda's SearchPoints. Check = Round[10^4 * obj] (2 -> well found).
_A2_D, _A2_S, _A2_C = 10, 2.0, 1.2345


def _a2_obj(x):
    r2 = np.sum((x - _A2_C) ** 2)
    return 1 - np.exp(-_A2_S * r2) + 1e-5 * np.sum(x ** 2)


_a2_bnds = [(-5.0, 5.0)] * _A2_D


def _a2_solve():
    return differential_evolution(_a2_obj, _a2_bnds, seed=1, tol=1e-9,
                                  popsize=40, maxiter=300, polish=True)


bench("A2 gaussian well (DE)", _a2_solve)
check("A2 gaussian well (DE)", int(np.floor(1e4 * _a2_solve().fun + 0.5)))


# ---- A3: modified-Ackley product-of-cosines (10-D) -- dual_annealing ----
# True global -1 at x=0 unreachable; both solvers plateau ~ -0.95..-0.98.
# Check = round(obj) = -1 for both.
from scipy.optimize import dual_annealing  # noqa: E402

_A3_D = 10


def _a3_obj(x):
    rms = np.sqrt(np.mean(x ** 2))
    return -np.exp(-0.2 * rms) * np.prod(np.cos(20 * x)) + 0.05 * np.sum(x ** 2)


_a3_bnds = [(-5.0, 5.0)] * _A3_D


def _a3_solve():
    return dual_annealing(_a3_obj, _a3_bnds, seed=1, maxiter=1000)


bench("A3 modified Ackley (SA)", _a3_solve)
check("A3 modified Ackley (SA)", int(np.floor(_a3_solve().fun + 0.5)))


# ---- A4: 15-D constrained minimax (Chebyshev) -- SLSQP on the epigraph ----
# Convex, unique optimum 0.125116. Matching character for NelderMead: a local
# constrained solver = scipy.optimize.minimize (SLSQP). vars z = [x_0..x_14, t].
from scipy.optimize import minimize  # noqa: E402

_A4_D = 15
_A4_S = np.sin(np.arange(1, _A4_D + 1))
_A4_IDX = np.arange(1, _A4_D + 1)


def _a4_cons():
    cons = [{'type': 'ineq', 'fun': lambda z: _A4_D - np.sum(z[:_A4_D] ** 2)},
            {'type': 'eq', 'fun': lambda z: np.sum(_A4_IDX * z[:_A4_D]) - 1}]
    for i in range(_A4_D):
        cons.append({'type': 'ineq', 'fun': (lambda z, i=i: z[-1] - (z[i] - _A4_S[i]))})
        cons.append({'type': 'ineq', 'fun': (lambda z, i=i: z[-1] - (_A4_S[i] - z[i]))})
    return cons


_a4_bnds = [(-2, 2)] * _A4_D + [(0, None)]
_a4_z0 = np.append(np.zeros(_A4_D), 1.0)


def _a4_solve():
    return minimize(lambda z: z[-1], _a4_z0, method='SLSQP', bounds=_a4_bnds,
                    constraints=_a4_cons(), options=dict(maxiter=500, ftol=1e-12))


bench("A4 minimax chebyshev (NM)", _a4_solve)
check("A4 minimax chebyshev (NM)", int(np.floor(1e4 * _a4_solve().x[-1] + 0.5)))


# ---- A5: Almgren-Chriss optimal liquidation (T=15) -- SLSQP single-start ----
# Convex; optimum 0.035206. Interior vars z = x[1..T-1]; x[0]=1, x[T]=0 pinned.
_A5_T, _A5_ETA, _A5_LAM, _A5_S2 = 15, 0.1, 0.05, 0.04


def _a5_expand(z):
    x = np.empty(_A5_T + 1); x[0] = 1.0; x[-1] = 0.0; x[1:_A5_T] = z
    return x


def _a5_obj(z):
    x = _a5_expand(z)
    return np.sum(_A5_ETA * np.abs(x[:-1] - x[1:]) ** 1.5 + _A5_LAM * _A5_S2 * x[:-1] ** 2)


_a5_cons = [{'type': 'ineq', 'fun': (lambda z, t=t: _a5_expand(z)[t] - _a5_expand(z)[t + 1])}
            for t in range(_A5_T)]
_a5_bnds = [(0, 1)] * (_A5_T - 1)
_a5_z0 = np.linspace(1, 0, _A5_T + 1)[1:_A5_T]


def _a5_solve():
    return minimize(_a5_obj, _a5_z0, method='SLSQP', bounds=_a5_bnds,
                    constraints=_a5_cons, options=dict(maxiter=500, ftol=1e-12))


bench("A5 liquidation (SA)", _a5_solve)
check("A5 liquidation (SA)", int(np.floor(1e6 * _a5_solve().fun + 0.5)))
