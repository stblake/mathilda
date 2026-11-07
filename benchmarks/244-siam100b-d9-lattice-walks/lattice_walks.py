"""Experiment 244 -- holdout B for SIAM Appendix D Problem 9 (Python column).

The ten transition matrices of lattice_walks.m: numpy.linalg.matrix_power at
machine precision (row of the start vertex of P^n, summed over the targets);
the 50-digit row pushes the start vector through n products v * P of mpmath
matrices at 60 digits, as experiment 184 (mpmath has no exact-rational or
sparse matrix type).  See lattice_walks.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
from fractions import Fraction
import itertools
import numpy as np
import mpmath

require(["numpy.linalg:matrix_power", "mpmath:matrix"])

def mkP(pts, nb):
    ix = {p: i for i, p in enumerate(pts)}
    P = np.full((len(pts), len(pts)), Fraction(0), dtype=object)
    for p in pts:
        l = nb(p)
        for q in l:
            P[ix[p], ix[q]] += Fraction(1, len(l))
    return P, ix
def grid_pts(d):
    return list(itertools.product(*[range(k) for k in d]))
def grid_nb(d, mv):
    return lambda p: [q for q in (tuple(a + b for a, b in zip(p, m)) for m in mv)
                      if all(0 <= t < k for t, k in zip(q, d))]
def tri_pts(n):
    return [(r, c) for r in range(n) for c in range(r + 1)]
def tri_nb(n):
    return lambda p: [(a, b) for (a, b) in ((p[0], p[1]-1), (p[0], p[1]+1), (p[0]-1, p[1]-1), (p[0]-1, p[1]),
                                            (p[0]+1, p[1]), (p[0]+1, p[1]+1)) if 0 <= b <= a < n]
def walk(pts, nb, s, tg, n):
    P, ix = mkP(pts, nb)
    return P, ix[s], [ix[t] for t in tg], n

king = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0)]
knight = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]
cube = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
def W(k):
    if k == 1: return walk(grid_pts([21]), grid_nb([21], [(1,), (-1,)]), (0,), [(i,) for i in range(11, 21)], 40)
    if k == 2: return walk(grid_pts([10, 10]), grid_nb([10, 10], [(1, 0), (-1, 0), (0, 1), (0, -1)]), (0, 0), [(9, j) for j in range(10)], 30)
    if k == 3: return walk(tri_pts(20), tri_nb(20), (0, 0), [(19, c) for c in range(20)], 40)
    if k == 4: return walk(grid_pts([6] * 3), grid_nb([6] * 3, cube), (0, 0, 0), [p for p in grid_pts([6] * 3) if p[0] == 5], 24)
    if k == 5: return walk(grid_pts([8, 8]), grid_nb([8, 8], knight), (0, 0), [(7, 7)], 20)
    if k == 6: return walk(grid_pts([8, 8]), grid_nb([8, 8], king), (0, 0), [(i, 7) for i in range(8)], 30)
    if k == 7: return walk(list(range(1024)), lambda v: [v ^ (1 << i) for i in range(10)], 0, [1023], 40)
    if k == 8: return walk(list(range(17)), lambda v: [(v + 1) % 17, (v - 1) % 17], 0, [0], 60)
    if k == 9:
        P = np.full((31, 31), Fraction(0), dtype=object)
        P[0, 0] = P[30, 30] = Fraction(1)
        for v in range(1, 30):
            P[v, v + 1] = Fraction(2, 5); P[v, v - 1] = Fraction(3, 5)
        return P, 15, [30], 200
    if k == 10: return walk(tri_pts(40), tri_nb(40), (0, 0), [(39, c) for c in range(40)], 100)

for k in range(1, 11):
    Pq, s, tg, n = W(k)
    Pf = Pq.astype(float)
    vrun("bd9.%d machine" % k, lambda: np.linalg.matrix_power(Pf, n)[s, tg].sum(), 10, 3, cap=60)
    def wp50():
        m = Pq.shape[0]
        P = mpmath.matrix(m, m)
        for i, j in zip(*np.nonzero(Pf)):
            q = Pq[i, j]
            P[int(i), int(j)] = mpmath.mpf(q.numerator) / q.denominator
        v = mpmath.matrix(1, m); v[s] = 1
        for _ in range(n):
            v = v * P
        return mpmath.fsum(v[0, j] for j in tg)
    vrun("bd9.%d wp50" % k, wp50, 48, cap=120, dps=60)
