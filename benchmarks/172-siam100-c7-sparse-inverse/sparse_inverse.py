"""Experiment 172 -- SIAM 100-Digit Challenge, Problem 7 (Python column).

A as a scipy.sparse matrix (diags), (A^-1)_11 = first component of
scipy.sparse.linalg.spsolve(A, e1) at machine precision -- at n = 2000.  At
n = 20000 spsolve (SuperLU) does not finish in 300 s with either the default
COLAMD or the symmetric MMD_AT_PLUS_A ordering (the 2^k bands fill the
factor in), and it cannot be interrupted by the per-case SIGALRM cap (C code),
so c7.machine uses the library's iterative solver for a symmetric positive
definite matrix instead: scipy.sparse.linalg.cg(A, e1, rtol=1e-14) with no
preconditioner argument.  mpmath has no
sparse solver: its only linear solver is the dense mpmath.lu_solve, so
  c7.wp50 (n = 20000)  SKIP -- a dense 20000^2 mpf matrix (400M objects) does
                       not fit in memory;
  c7.wp50 n=2000       mpmath.lu_solve on the dense 2000 x 2000 matrix at
                       dps 60 (O(n^3) mpf operations; expected to hit the cap).
The matrix is built outside the timed call.  See sparse_inverse.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require, mp
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import mpmath

require(["scipy.sparse.linalg:spsolve", "scipy.sparse.linalg:cg", "mpmath:lu_solve"])


def primes(n):
    N = 300000
    s = np.ones(N, bool); s[:2] = False
    for i in range(2, int(N ** .5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0][:n].astype(float)


def mat(n):
    offs = [2 ** k for k in range(n.bit_length()) if 2 ** k < n]
    diags = [primes(n)] + [np.ones(n - o) for o in offs] * 2
    pos = [0] + offs + [-o for o in offs]
    return sp.diags(diags, pos, format="csc")


A20, A2 = mat(20000), mat(2000)
e20 = np.zeros(20000); e20[0] = 1
e2 = np.zeros(2000); e2[0] = 1

vrun("c7.machine", lambda: sla.cg(A20, e20, rtol=1e-14, atol=0.0, maxiter=20000)[0][0], 10)
skip("c7.wp50", "mpmath.lu_solve (dense 20000^2)")
vrun("c7.machine n=2000", lambda: sla.spsolve(A2, e2)[0], 10, 3)

mp.dps = 60
A2w = mpmath.matrix(A2.toarray().tolist())
e2w = mpmath.matrix([1] + [0] * 1999)
vrun("c7.wp50 n=2000", lambda: mpmath.lu_solve(A2w, e2w)[0], 48, dps=60)
