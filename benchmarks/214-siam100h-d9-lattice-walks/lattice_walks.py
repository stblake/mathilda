"""Experiment 214 -- SIAM 100-digit HOLDOUT A, family D9 (Python column).

The ten walks of lattice_walks.m.  machine: numpy.linalg.matrix_power of the
float transition matrix, start row summed over the targets; wp50: the start
vector pushed through `steps` products v * P of mpmath matrices at 60 digits
(mpmath has no exact-rational matrix type), as in 184.  The matrices are
built (as exact fractions, then converted) before the timed cases."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
from fractions import Fraction
import numpy as np
import mpmath

require(["numpy.linalg:matrix_power", "mpmath:matrix"])

def tmat(pts, nb, mode="plain"):
    ix = {p: i for i, p in enumerate(pts)}
    n = len(pts)
    P = np.full((n, n), Fraction(0), dtype=object)
    for p in pts:
        l, i = nb(p), ix[p]
        if mode == "lazy":
            P[i, i] = Fraction(1, 2)
            for q in l: P[i, ix[q]] += Fraction(1, 2 * len(l))
        elif mode != "plain" and mode(p):
            P[i, i] = Fraction(1)
        else:
            for q in l: P[i, ix[q]] += Fraction(1, len(l))
    return P
def pos(pts, test): return [i for i, p in enumerate(pts) if test(p)]
tri_pts = lambda n: [(r, c) for r in range(n) for c in range(r + 1)]
def tri_nb(n):
    return lambda p: [q for q in ((p[0], p[1]-1), (p[0], p[1]+1), (p[0]-1, p[1]-1), (p[0]-1, p[1]),
                                  (p[0]+1, p[1]), (p[0]+1, p[1]+1)) if 0 <= q[1] <= q[0] < n]
grid_pts = lambda n: [(i, j) for i in range(1, n+1) for j in range(1, n+1)]
rook = [(1, 0), (-1, 0), (0, 1), (0, -1)]
king = rook + [(1, 1), (1, -1), (-1, 1), (-1, -1)]
def grid_nb(n, ds):
    return lambda p: [(p[0]+a, p[1]+b) for a, b in ds if 1 <= p[0]+a <= n and 1 <= p[1]+b <= n]
cube_pts = lambda n: [(i, j, k) for i in range(1, n+1) for j in range(1, n+1) for k in range(1, n+1)]
def cube_nb(n):
    ds = [(1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)]
    return lambda p: [q for q in (tuple(a+b for a, b in zip(p, d)) for d in ds) if min(q) >= 1 and max(q) <= n]
hyp_nb = lambda d: (lambda p: [p ^ (1 << i) for i in range(d)])

T10, T20, T30, T40 = tri_pts(10), tri_pts(20), tri_pts(30), tri_pts(40)
G20, G12, C6, H8 = grid_pts(20), grid_pts(12), cube_pts(6), list(range(256))
prob = {
  1: (tmat(T10, tri_nb(10)), 0, 20, pos(T10, lambda p: p[0] == 9)),
  2: (tmat(T30, tri_nb(30)), 0, 60, [0]),
  3: (tmat(T20, tri_nb(20)), pos(T20, lambda p: p == (12, 6))[0], 50,
      pos(T20, lambda p: p[1] == 0 or p[1] == p[0] or p[0] == 19)),
  4: (tmat(G20, grid_nb(20, rook)), 0, 60, pos(G20, lambda p: p[0] == 20)),
  5: (tmat(G12, grid_nb(12, king)), 0, 30, pos(G12, lambda p: p[0] == p[1])),
  6: (tmat(C6, cube_nb(6)), 0, 41, [len(C6) - 1]),
  7: (tmat(H8, hyp_nb(8)), 0, 30, [255]),
  8: (tmat(T30, tri_nb(30), "lazy"), 0, 60, pos(T30, lambda p: p[0] == 29)),
  9: (tmat(T30, tri_nb(30), lambda p: p[0] == 29), 0, 60, pos(T30, lambda p: p[0] == 29)),
  10: (tmat(T40, tri_nb(40)), 0, 80, pos(T40, lambda p: p[0] == 39)),
}
PF = {k: pr[0].astype(float) for k, pr in prob.items()}
def walk_m(k):
    _, s, n, tg = prob[k]
    return np.linalg.matrix_power(PF[k], n)[s, tg].sum()
def walk_50(pr):
    Pq, s, k, tg = pr
    n = Pq.shape[0]
    P = mpmath.matrix(n, n)
    for i, j in zip(*np.nonzero(Pq != 0)):
        q = Pq[i, j]
        P[int(i), int(j)] = mpmath.mpf(q.numerator) / q.denominator
    v = mpmath.matrix(1, n); v[s] = 1
    for _ in range(k):
        v = v * P
    return mpmath.fsum(v[0, j] for j in tg)
for k in range(1, 11):
    vrun("hd9.%d machine" % k, lambda: walk_m(k), 10, 3, cap=30)
    vrun("hd9.%d wp50" % k, lambda: walk_50(prob[k]), 48, cap=90, dps=60)
