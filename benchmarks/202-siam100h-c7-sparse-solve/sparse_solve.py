"""Experiment 202 -- SIAM 100-Digit Challenge, holdout A for Problem 7 (Python column).

The same ten integer sparse matrices as sparse_solve.m, as scipy.sparse
(diags); machine rows = scipy.sparse.linalg.spsolve (SuperLU) of the
matrix against the right-hand side, then the requested entry / total.
mpmath has no sparse solver: its only linear solver is the dense
mpmath.lu_solve, so the wp50 rows are
  n = 1000 (hc7.5, hc7.7, hc7.9)  mpmath.lu_solve on the dense matrix at dps 60
                                  (O(n^3) mpf operations; capped at 60 s);
  n >= 3000 (the rest)            SKIP -- a dense n^2 mpf matrix.
The matrices are built outside the timed call.  See sparse_solve.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require, mp
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import mpmath

require(["scipy.sparse.linalg:spsolve", "mpmath:lu_solve"])


def primes(n):
    N = 250000
    s = np.ones(N, bool); s[:2] = False
    for i in range(2, int(N ** .5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0][:n].astype(float)


def mat(n, offs, d):
    diags = [np.asarray(d, float)] + [np.ones(n - o) for o in offs] * 2
    return sp.diags(diags, [0] + list(offs) + [-o for o in offs], format="csc")


def pow2(n):
    return [2 ** k for k in range(n.bit_length()) if 2 ** k < n]


def fib(n):
    f, a, b = [], 1, 2
    while a < n:
        f.append(a); a, b = b, a + b
    return f


def lap2(m):
    n = m * m
    r = np.array([0.0 if (i % m) == 0 else -1.0 for i in range(1, n)])
    return sp.diags([4.0 * np.ones(n), r, r, -np.ones(n - m), -np.ones(n - m)], [0, 1, -1, m, -m], format="csc")


def lap3(m):
    n = m ** 3
    r1 = np.array([0.0 if i % m == 0 else -1.0 for i in range(1, n)])
    rm = np.array([0.0 if (i % (m * m) > m * m - m or i % (m * m) == 0) else -1.0 for i in range(1, n - m + 1)])
    return sp.diags([6.0 * np.ones(n), r1, r1, rm, rm, -np.ones(n - m * m), -np.ones(n - m * m)],
                    [0, 1, -1, m, -m, m * m, -m * m], format="csc")


P = primes(20000)
M = {1: mat(20000, [1], P), 2: mat(10000, [1, 2], P[:10000]), 3: mat(5000, [3 ** k for k in range(8)], P[:5000]),
     4: mat(3000, fib(3000), P[:3000]), 5: mat(1000, pow2(1000), np.arange(1, 1001)),
     6: mat(3000, pow2(3000), P[:3000]), 7: mat(1000, pow2(1000), P[:1000]), 8: lap2(60),
     9: mat(1000, pow2(1000), P[:1000]), 10: lap3(20)}


def unit(n, i):
    b = np.zeros(n); b[i] = 1.0
    return b


# (rhs, functional) per case
F = {1: ("e1", lambda x: x[0]), 2: ("e1", lambda x: x[0]), 3: ("e1", lambda x: x[0]), 4: ("e1", lambda x: x[0]),
     5: ("e1", lambda x: x[0]), 6: ("en", lambda x: x[-1]), 7: ("e1", lambda x: x[1]), 8: ("e1", lambda x: x[0]),
     9: ("ones", lambda x: sum(x)), 10: ("e1", lambda x: x[0])}


def rhs(n, kind):
    return unit(n, 0) if kind == "e1" else unit(n, n - 1) if kind == "en" else np.ones(n)


for k in range(1, 11):
    A = M[k]; n = A.shape[0]; kind, fn = F[k]
    b = rhs(n, kind)
    vrun("hc7.%d machine" % k, lambda A=A, b=b, fn=fn: float(fn(sla.spsolve(A, b))), 10, 3 if n <= 10000 else 1, cap=60)
    if n <= 1000:
        mp.dps = 60
        Aw = mpmath.matrix(A.toarray().astype(int).tolist())
        bw = mpmath.matrix([int(v) for v in b])
        vrun("hc7.%d wp50" % k, lambda Aw=Aw, bw=bw, fn=fn: fn(list(mpmath.lu_solve(Aw, bw))), 48, cap=60, dps=60)
    else:
        skip("hc7.%d wp50" % k, "mpmath.lu_solve (dense %d^2)" % n)
