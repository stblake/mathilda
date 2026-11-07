"""Experiment 224 -- SIAM 100-Digit Challenge holdout A, family hd19 (Python column).

Finite variable-thickness design problems (see design_variants.m): minimise
the compliance J(c) = F.U, K(c) U = F, over lo <= c_e <= hi with
sum(c) = mean * ncells.  scipy.optimize.minimize (SLSQP, bounds + equality
constraint, numerical gradient) at machine precision, as in exp 194; mpmath
has no constrained optimiser, so every wp50 row is a SKIP."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import numpy as np
from scipy.optimize import minimize

require(["scipy.optimize:minimize", "numpy.linalg:solve"])

KQ = np.array([[4, -1, -2, -1], [-1, 4, -1, -2], [-2, -1, 4, -1], [-1, -2, -1, 4]]) / 6.0
# (dim, nx, ny, lo, hi, mean, load f(x))
PROBS = [
    (2, 4, 4, 1, 10, 11 / 2, lambda x: 1.0), (2, 4, 4, 1, 100, 103 / 4, lambda x: 1.0),
    (2, 6, 6, 1, 100, 101 / 2, lambda x: 1.0), (2, 5, 5, 1, 100, 101 / 2, lambda x: 1.0),
    (1, 8, 0, 1, 100, 101 / 2, lambda x: 1.0), (1, 16, 0, 1, 4, 5 / 2, lambda x: 1.0),
    (2, 4, 4, 1, 100, 101 / 2, lambda x: x), (2, 8, 4, 1, 100, 101 / 2, lambda x: 1.0),
    (2, 4, 4, 1, 2, 3 / 2, lambda x: 1.0), (2, 6, 6, 1, 1000, 1003 / 4, lambda x: 1.0)]


def mesh(dim, nx, ny, f):
    if dim == 1:
        h = 1.0 / nx
        node = lambda i: i - 1 if 1 <= i <= nx - 1 else -1
        els = [[node(a), node(a + 1)] for a in range(nx)]
        return els, np.array([[1, -1], [-1, 1]]) / h, np.array([h * f(h * i) for i in range(1, nx)]), nx - 1
    h = 1.0 / ny
    node = lambda i, j: (i - 1) * (ny - 1) + (j - 1) if 1 <= i <= nx - 1 and 1 <= j <= ny - 1 else -1
    els = [[node(a, b), node(a + 1, b), node(a + 1, b + 1), node(a, b + 1)] for a in range(nx) for b in range(ny)]
    F = np.array([h * h * f(h * i) for i in range(1, nx) for j in range(1, ny)])
    return els, KQ, F, (nx - 1) * (ny - 1)


def compliance(ms, c):
    els, loc, F, m = ms
    K = np.zeros((m, m))
    for e, idx in enumerate(els):
        for r, ir in enumerate(idx):
            if ir < 0:
                continue
            for s, js in enumerate(idx):
                if js >= 0:
                    K[ir, js] += c[e] * loc[r, s]
    return F @ np.linalg.solve(K, F)


def design_min(p):
    dim, nx, ny, lo, hi, mean, f = p
    ms = mesh(dim, nx, ny, f)
    ne = len(ms[0])
    res = minimize(lambda c: compliance(ms, c), np.full(ne, mean), method="SLSQP",
                   bounds=[(lo, hi)] * ne,
                   constraints=[{"type": "eq", "fun": lambda c: c.sum() - mean * ne}],
                   options={"ftol": 1e-15, "maxiter": 2000})
    return res.fun


for k, p in enumerate(PROBS, 1):
    vrun("hd19.%d machine" % k, lambda: design_min(p), 10, cap=90)
    skip("hd19.%d wp50" % k, "ConstrainedMinimize(mpmath)")
