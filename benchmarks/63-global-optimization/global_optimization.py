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
