"""Experiment 232 -- SIAM 100-digit HOLDOUT B for Problem 7 (Python column).

The same ten matrices as sparse_solves.m, built as scipy.sparse matrices
outside the timed call; the case is scipy.sparse.linalg.spsolve(A, b)[j]
(SuperLU, default COLAMD ordering).  Exception, as in experiment 172: for
bc7.3 (n = 20000, the 2^k-like Fibonacci bands fill the factor in) spsolve
does not finish in the cap and cannot be interrupted (C code), so that row
uses the library's solver for symmetric positive definite matrices,
scipy.sparse.linalg.cg(A, b, rtol=1e-14) with no preconditioner argument.
mpmath has no sparse solver: its only linear solver is the dense
mpmath.lu_solve, so the wp50 rows are mpmath.lu_solve on the dense matrix at
dps 60 for n <= 2047 (bc7.1, 4, 5, 7, 10; expected to hit the cap) and SKIP
for the large n (a dense n^2 mpf matrix does not fit).  See sparse_solves.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require, mp
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import mpmath

require(["scipy.sparse.linalg:spsolve", "scipy.sparse.linalg:cg", "mpmath:lu_solve"])


def primes(n):
    N = 20 * n + 100
    s = np.ones(N, bool); s[:2] = False
    for i in range(2, int(N ** .5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0][:n].astype(float)


def prime_bands(n, offs, sup=1.0, sub=1.0):
    d = [primes(n)] + [np.full(n - o, sup) for o in offs] + [np.full(n - o, sub) for o in offs]
    return sp.diags(d, [0] + list(offs) + [-o for o in offs], format="csc")


def lap1(m):
    return sp.diags([np.full(m, 2.0), np.full(m - 1, -1.0), np.full(m - 1, -1.0)], [0, 1, -1])


def lap2(m):
    I = sp.identity(m)
    return sp.kron(lap1(m), I) + sp.kron(I, lap1(m))


def lap3(m):
    return sp.kron(lap1(m), sp.identity(m * m)) + sp.kron(sp.identity(m), lap2(m))


def cayley(p):
    r, c = [], []
    for v in range(1, p):
        for g in (2, 3, 5):
            w = g * v % p
            r += [v, w]; c += [w, v]
    return (sp.csc_matrix((np.ones(len(r)), (r, c)), shape=(p, p)) + 7 * sp.identity(p)).tocsc()


def tree(n):
    r, c, v = [], [], []
    for i in range(1, n + 1):
        r.append(i - 1); c.append(i - 1)
        v.append(3.0 if i == 1 else (4.0 if 2 * i <= n else 2.0))
    for i in range(1, (n - 1) // 2 + 1):
        for ch in (2 * i, 2 * i + 1):
            r += [i - 1, ch - 1]; c += [ch - 1, i - 1]; v += [-1.0, -1.0]
    return sp.csc_matrix((v, (r, c)), shape=(n, n))


fib = [1, 2]
while fib[-1] + fib[-2] < 20000:
    fib.append(fib[-1] + fib[-2])
A = {1: prime_bands(1000, [1]),
     2: (lap2(100) + sp.identity(10 ** 4)).tocsc(),
     3: prime_bands(20000, fib),
     4: prime_bands(2000, [j * j for j in range(1, 45)]),
     5: prime_bands(2000, [2 ** j for j in range(11)], 1.0, -1.0),
     6: prime_bands(5000, [2 ** j for j in range(13)]),
     7: cayley(2003),
     8: lap3(30).tocsc(),
     9: (lap2(150) - 0.5 * sp.identity(22500)).tocsc(),
     10: tree(2047)}
IDX = {1: 1, 2: 4950, 3: 1, 4: 1, 5: 1, 6: 1, 7: 2, 8: 1, 9: 11175, 10: 1}   # 1-based
DENSE_OK = (1, 4, 5, 7, 10)


def rhs(k):
    n = A[k].shape[0]
    if k == 6:
        return np.ones(n)
    b = np.zeros(n); b[IDX[k] - 1] = 1.0
    return b


for k in range(1, 11):
    j = IDX[k] - 1
    Ak, bk = A[k], rhs(k)
    if k == 3:
        vrun("bc7.3 machine", lambda: sla.cg(Ak, bk, rtol=1e-14, atol=0.0, maxiter=20000)[0][j], 10, cap=90)
    else:
        vrun("bc7.%d machine" % k, lambda: sla.spsolve(Ak, bk)[j], 10, cap=90)
    if k in DENSE_OK:
        mp.dps = 60
        Aw = mpmath.matrix(Ak.toarray().tolist())
        bw = mpmath.matrix(bk.tolist())
        vrun("bc7.%d wp50" % k, lambda: mpmath.lu_solve(Aw, bw)[j], 48, cap=75, dps=60)
        del Aw
    else:
        skip("bc7.%d wp50" % k, "mpmath.lu_solve (dense %d^2)" % Ak.shape[0])
