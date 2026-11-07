"""Experiment 184 -- SIAM 100-Digit Challenge, Appendix D Problem 9 (Python column).

The 465 x 465 transition matrix P of the walk: numpy.linalg.matrix_power at
machine precision (first row of P^60, summed over the bottom row); the 50-digit
row pushes the start vector through 60 products v * P of mpmath matrices at
60 digits (mpmath has no exact-rational matrix type; an object-dtype numpy
matrix of fractions.Fraction does the same 60 exact products in ~95 s).
See triangular_walk.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require, mp
from fractions import Fraction
import numpy as np
import mpmath

require(["numpy.linalg:matrix_power", "mpmath:matrix"])

nT = 30
idx = lambda r, c: r * (r + 1) // 2 + c
def nbrs(r, c):
    return [(a, b) for (a, b) in ((r, c-1), (r, c+1), (r-1, c-1), (r-1, c), (r+1, c), (r+1, c+1))
            if 0 <= b <= a < nT]
n = nT * (nT + 1) // 2
Pq = np.full((n, n), Fraction(0), dtype=object)
for r in range(nT):
    for c in range(r + 1):
        l = nbrs(r, c)
        for (a, b) in l:
            Pq[idx(r, c), idx(a, b)] = Fraction(1, len(l))
Pf = Pq.astype(float)
bottom = [idx(nT - 1, c) for c in range(nT)]

vrun("d9.machine", lambda: np.linalg.matrix_power(Pf, 60)[0, bottom].sum(), 10, 3)

def wp50():
    P = mpmath.matrix(n, n)
    for i, j in zip(*np.nonzero(Pf)):
        P[int(i), int(j)] = mpmath.mpf(1) / round(1 / Pf[i, j])
    v = mpmath.matrix(1, n); v[0] = 1
    for _ in range(60):
        v = v * P
    return mpmath.fsum(v[0, j] for j in bottom)
vrun("d9.wp50", wp50, 48, dps=60)
