"""Experiment 194 -- SIAM 100-Digit Challenge, App. D Problem 19 (Python column).

The finite proxy of the relaxed optimal-design problem (see optimal_design.m):
n = 8 (and n = 4 for the reduced row) cells per side, Q1 elements, minimise the compliance J_8(c) over
1 <= c_e <= 100 with sum(c) = 50.5 n^2.  scipy.optimize.minimize (SLSQP,
bounds + equality constraint, numerical gradient) at machine precision;
mpmath has no constrained optimiser, so the wp50 row is a SKIP."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import numpy as np
from scipy.optimize import minimize

require(["scipy.optimize:minimize", "numpy.linalg:solve"])

KQ = np.array([[4, -1, -2, -1], [-1, 4, -1, -2], [-2, -1, 4, -1], [-1, -2, -1, 4]]) / 6.0


def compliance(n, c):
    node = lambda i, j: (i - 1) * (n - 1) + (j - 1) if 1 <= i <= n - 1 and 1 <= j <= n - 1 else -1
    N = (n - 1) ** 2
    K = np.zeros((N, N))
    for a in range(n):
        for b in range(n):
            idx = [node(a, b), node(a + 1, b), node(a + 1, b + 1), node(a, b + 1)]
            for r in range(4):
                for s in range(4):
                    if idx[r] >= 0 and idx[s] >= 0:
                        K[idx[r], idx[s]] += c[a * n + b] * KQ[r, s]
    F = np.full(N, 1.0 / n ** 2)
    return F @ np.linalg.solve(K, F)


def design_min(n):
    res = minimize(lambda c: compliance(n, c), np.full(n * n, 50.5), method="SLSQP",
                   bounds=[(1, 100)] * (n * n),
                   constraints=[{"type": "eq", "fun": lambda c: c.sum() - 50.5 * n * n}],
                   options={"ftol": 1e-15, "maxiter": 1000})
    return res.fun


vrun("d19.machine", lambda: design_min(8), 10)
vrun("d19.machine-reduced", lambda: design_min(4), 10)
skip("d19.wp50", "ConstrainedMinimize(mpmath)")
