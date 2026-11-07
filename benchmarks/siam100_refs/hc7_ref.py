"""Holdout A, family C7 (experiment 202): references for hc7.1 .. hc7.10.

Method A (this script): mixed-precision iterative refinement.  A is an
INTEGER sparse matrix; the correction is solved in doubles with a sparse LU
(scipy splu, which also handles the indefinite hc7.5), and the residual
b - A x is computed EXACTLY in Python integers with x held as an integer scaled
by 2^PB, so every sweep gains ~10-15 digits until x is exact to ~2^-PB.
Method B: siam100_refs/hc7_ref.wl (wolframscript LinearSolve on the same
SparseArray at WorkingPrecision 80 -- a different algorithm and arithmetic).
Method C (independent closed forms, printed here too):
  hc7.1  the continued fraction 1/(a1 - 1/(a2 - 1/(a3 - ...))) of a tridiagonal
         matrix with unit off-diagonals, evaluated backwards at dps 110;
  hc7.8  the eigen-sum (2/(m+1))^2 sum_{j,k} sin^2(j pi/(m+1)) sin^2(k pi/(m+1))
         / (4 - 2 cos(j pi/(m+1)) - 2 cos(k pi/(m+1)));
  hc7.10 the 3-D eigen-sum likewise.
Usage: python3.11 hc7_ref.py
"""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from mpmath import mp, mpf, nstr, sin, cos, pi, fsum

mp.dps = 110
PB = 360
S = 1 << PB


def primes(n):
    N = 250000
    s = np.ones(N, bool); s[:2] = False
    for i in range(2, int(N ** .5) + 1):
        if s[i]:
            s[i * i::i] = False
    return [int(p) for p in np.nonzero(s)[0][:n]]


def band_matrix(n, offs, diag):
    ent = {}
    for i in range(n):
        ent[(i, i)] = diag[i]
    for o in offs:
        for i in range(n - o):
            ent[(i, i + o)] = 1
            ent[(i + o, i)] = 1
    return n, ent


def lap(m, dim):
    n = m ** dim
    ent = {}
    for i in range(n):
        ent[(i, i)] = 2 * dim
        c = [(i // m ** d) % m for d in range(dim)]
        for d in range(dim):
            if c[d] < m - 1:
                j = i + m ** d
                ent[(i, j)] = -1
                ent[(j, i)] = -1
    return n, ent


def pow2(n):
    o, k = [], 1
    while k < n:
        o.append(k); k *= 2
    return o


def fib(n):
    f, a, b = [], 1, 2
    while a < n:
        f.append(a); a, b = b, a + b
    return f


def refine(n, ent, b):
    rows = [[] for _ in range(n)]
    for (i, j), v in ent.items():
        rows[i].append((j, v))
    I = [k[0] for k in ent]; J = [k[1] for k in ent]; V = [float(v) for v in ent.values()]
    lu = sla.splu(sp.csc_matrix((V, (I, J)), shape=(n, n)))
    X = [0] * n
    Bs = [bi * S for bi in b]
    for it in range(40):
        R = [Bs[i] - sum(v * X[j] for j, v in rows[i]) for i in range(n)]
        m = max(abs(r) for r in R)
        if m == 0:
            break
        e = max(0, m.bit_length() - 60)
        d = lu.solve(np.array([float(r >> e) if e else float(r) for r in R]))
        X = [X[i] + (int(round(d[i] * 2.0 ** 40)) << e >> 40) for i in range(n)]
        if m.bit_length() < 40:
            break
    print("  sweeps", it, "final residual bits", m.bit_length(), file=sys.stderr)
    return [mpf(x) / S for x in X]


cases = {}
P20 = primes(20000)
cases[1] = (band_matrix(20000, [1], P20), "e1", lambda x: x[0])
cases[2] = (band_matrix(10000, [1, 2], P20[:10000]), "e1", lambda x: x[0])
cases[3] = (band_matrix(5000, [3 ** k for k in range(8)], P20[:5000]), "e1", lambda x: x[0])
cases[4] = (band_matrix(3000, fib(3000), P20[:3000]), "e1", lambda x: x[0])
cases[5] = (band_matrix(1000, pow2(1000), list(range(1, 1001))), "e1", lambda x: x[0])
cases[6] = (band_matrix(3000, pow2(3000), P20[:3000]), "en", lambda x: x[-1])
cases[7] = (band_matrix(1000, pow2(1000), P20[:1000]), "e1", lambda x: x[1])
cases[8] = (lap(60, 2), "e1", lambda x: x[0])
cases[9] = (band_matrix(1000, pow2(1000), P20[:1000]), "ones", lambda x: fsum(x))
cases[10] = (lap(20, 3), "e1", lambda x: x[0])

for k in range(1, 11):
    (n, ent), bk, fn = cases[k]
    b = [0] * n
    if bk == "e1": b[0] = 1
    elif bk == "en": b[-1] = 1
    else: b = [1] * n
    print("hc7.%d" % k, nstr(fn(refine(n, ent, b)), 75), flush=True)

# closed forms
a = P20
v = mpf(a[-1])
for i in range(len(a) - 2, -1, -1):
    v = a[i] - 1 / v
print("C hc7.1 continued fraction", nstr(1 / v, 75))
m = 60; h = pi / (m + 1)
s2 = [sin(j * h) ** 2 for j in range(1, m + 1)]; cs = [cos(j * h) for j in range(1, m + 1)]
print("C hc7.8 eigen-sum", nstr((2 / mpf(m + 1)) ** 2 * fsum(s2[j] * s2[k] / (4 - 2 * cs[j] - 2 * cs[k])
                                                         for j in range(m) for k in range(m)), 75))
m = 20; h = pi / (m + 1)
s2 = [sin(j * h) ** 2 for j in range(1, m + 1)]; cs = [cos(j * h) for j in range(1, m + 1)]
print("C hc7.10 eigen-sum", nstr((2 / mpf(m + 1)) ** 3 * fsum(s2[i] * s2[j] * s2[k] / (6 - 2 * cs[i] - 2 * cs[j] - 2 * cs[k])
                                                          for i in range(m) for j in range(m) for k in range(m)), 75))
