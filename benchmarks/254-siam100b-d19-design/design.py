"""Experiment 254 -- SIAM 100-digit SEALED HOLDOUT B, family D19 (Python column).

The ten finite convex design problems of design.m (the same element lists,
built here by the same rules).  scipy.optimize.minimize (SLSQP, bounds + the
equality constraint, numerical gradient: the objective is a black box, as in
exp 194) at machine precision; mpmath has no constrained optimiser, so every
wp50 row is a SKIP."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import minimize

require(["scipy.optimize:minimize", "numpy.linalg:solve"])

KQ = [[4, -1, -2, -1], [-1, 4, -1, -2], [-2, -1, 4, -1], [-1, -2, -1, 4]]   # x 1/6, corners ccw


def bar(n, clamp_right):
    h = Fr(1, n)
    dof = {i: i - 1 for i in range(1, n + (0 if clamp_right else 1))}
    N = len(dof)
    F = [Fr(0)] * N
    els = []
    for e in range(n):
        nodes = [e, e + 1]
        els.append(([dof.get(v, -1) for v in nodes], [[1 / h, -1 / h], [-1 / h, 1 / h]]))
        for v in nodes:
            if v in dof:
                F[dof[v]] += h / 2
    return N, els, F


def q1(nx, ny, h, cell_ok, node_free, point=None):
    cells = [(a, b) for a in range(nx) for b in range(ny) if cell_ok(a, b)]
    used = sorted({v for a, b in cells for v in ((a, b), (a + 1, b), (a + 1, b + 1), (a, b + 1))})
    free = [v for v in used if node_free(*v)]
    dof = {v: i for i, v in enumerate(free)}
    F = [Fr(0)] * len(free)
    els = []
    for a, b in cells:
        nodes = [(a, b), (a + 1, b), (a + 1, b + 1), (a, b + 1)]
        els.append(([dof.get(v, -1) for v in nodes], [[Fr(k, 6) for k in r] for r in KQ]))
        if point is None:
            for v in nodes:
                if v in dof:
                    F[dof[v]] += h * h / 4
    if point is not None:
        F[dof[point]] = Fr(1)
    return len(free), els, F


def interior(nx, ny):
    return lambda i, j: 0 < i < nx and 0 < j < ny


def p1cell(n):
    """each cell = two P1 right triangles split along its (a,b)-(a+1,b+1) diagonal, ONE c per
    cell: the summed cell matrix (1/2)[[2,-1,0,-1],[-1,2,-1,0],[0,-1,2,-1],[-1,0,-1,2]] (the
    5-point stencil) and the load h^2 (1/3, 1/6, 1/3, 1/6) on the corners (a,b),(a+1,b),
    (a+1,b+1),(a,b+1)"""
    h = Fr(1, n)
    dof = {(i, j): (i - 1) * (n - 1) + (j - 1) for i in range(1, n) for j in range(1, n)}
    F = [Fr(0)] * len(dof)
    Kc = [[Fr(x, 2) for x in r] for r in ([2, -1, 0, -1], [-1, 2, -1, 0], [0, -1, 2, -1], [-1, 0, -1, 2])]
    w = [h * h / 3, h * h / 6, h * h / 3, h * h / 6]
    els = []
    for a in range(n):
        for b in range(n):
            nodes = [(a, b), (a + 1, b), (a + 1, b + 1), (a, b + 1)]
            els.append(([dof.get(v, -1) for v in nodes], Kc))
            for v, wv in zip(nodes, w):
                if v in dof:
                    F[dof[v]] += wv
    return len(dof), els, F


def grid_graph(m):
    nodes = [(i, j) for i in range(m) for j in range(m)]
    dof = {v: k - 1 for k, v in enumerate(nodes) if k > 0}       # node (0,0) grounded
    els = []
    for i in range(m):
        for j in range(m):
            for nb in ((i + 1, j), (i, j + 1)):
                if nb[0] < m and nb[1] < m:
                    els.append(([dof.get((i, j), -1), dof.get(nb, -1)], [[Fr(1), Fr(-1)], [Fr(-1), Fr(1)]]))
    F = [Fr(0)] * len(dof)
    F[dof[(1, 2)]] = Fr(1)                                       # unit current in at (1,2)
    return len(dof), els, F


def problem(k):
    if k == 1: return bar(8, False) + (1, 10, 32)
    if k == 2: return bar(12, True) + (10, 100, 606)
    if k == 3: return q1(4, 4, Fr(1, 4), lambda a, b: True, interior(4, 4)) + (1, 10, 48)
    if k == 4: return q1(8, 4, Fr(1, 4), lambda a, b: True, interior(8, 4)) + (1, 100, 1616)
    if k == 5: return q1(6, 6, Fr(1, 6), lambda a, b: True, interior(6, 6), point=(3, 3)) + (1, 100, 1818)
    if k == 6: return q1(6, 6, Fr(1, 6), lambda a, b: True, lambda i, j: not (j == 0 and i <= 2)) + (1, 100, 360)
    if k == 7: return grid_graph(4) + (1, 10, 60)
    if k == 8: return p1cell(4) + (1, 100, 808)
    if k == 9:
        ok = lambda a, b: not (a >= 3 and b >= 3)
        cell = lambda i, j: all(ok(a, b) and 0 <= a < 6 and 0 <= b < 6 for a in (i - 1, i) for b in (j - 1, j))
        return q1(6, 6, Fr(1, 6), ok, cell) + (1, 100, 1350)
    if k == 10: return q1(10, 10, Fr(1, 10), lambda a, b: True, interior(10, 10)) + (1, 100, 5050)



def design_min(k):
    N, els, F, lo, hi, V = problem(k)
    els = [(np.array(idx), np.array([[float(x) for x in r] for r in Ke])) for idx, Ke in els]
    F = np.array([float(x) for x in F])
    m = len(els)

    def J(c):
        K = np.zeros((N, N))
        for (idx, Ke), ce in zip(els, c):
            ok = idx >= 0
            K[np.ix_(idx[ok], idx[ok])] += ce * Ke[np.ix_(ok, ok)]
        return F @ np.linalg.solve(K, F)

    res = minimize(J, np.full(m, V / m), method="SLSQP", bounds=[(lo, hi)] * m,
                   constraints=[{"type": "eq", "fun": lambda c: c.sum() - V}],
                   options={"ftol": 1e-15, "maxiter": 1000})
    return res.fun


for k in range(1, 11):
    vrun("bd19.%d machine" % k, lambda: design_min(k), 10, cap=90)
    skip("bd19.%d wp50" % k, "ConstrainedMinimize(mpmath)")
