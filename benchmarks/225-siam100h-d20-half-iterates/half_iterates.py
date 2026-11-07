"""Experiment 225 -- SIAM 100-Digit Challenge holdout A, family hd20 (Python column).

Regular half-iterates by the truncated Schroeder system (see half_iterates.m;
equation j divided by lam^j, as there).
Taylor coefficients of f by mpmath.taylor (the only Taylor tool in the
stack), powers f^k by truncated convolution (glue), then
  machine: numpy.linalg.solve + scipy.optimize.brentq on [a, b],
  wp50:    mpmath.lu_solve + mpmath.findroot (dps 60).
Bracket [a, b] = [z0, 4 z0] (hd20.9: [-1/2, -1/8]; hd20.10: [0.4, 0.7]), the
same as the Maxima / FriCAS halves.  hd20.10 is exp 195's Robbins system for
base 2 at n = 16."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require, mp
import numpy as np
import mpmath
from scipy.optimize import brentq

require(["mpmath:taylor", "numpy.linalg:solve", "scipy.optimize:brentq", "mpmath:lu_solve", "mpmath:findroot"])

S2 = lambda: mpmath.sqrt(2)
# (f, lam, z0, n machine, n wp50, bracket, post)
PROBS = [
    (lambda z: 16 * z / (4 + 3 * z), 4, (1, 2), 25, 90, None, lambda v: v),
    (lambda z: (4 * z - 4 * z ** 3) / (1 - 6 * z ** 2 + z ** 4), 4, (1, 10), 25, 90, None, lambda v: v),
    (lambda z: 4 * z * (1 - z), 4, (1, 10), 25, 90, None, lambda v: v),
    (lambda z: 4 * z * (1 + z), 4, (1, 10), 25, 90, None, lambda v: v),
    (lambda z: 1 - (4 * (1 - z) ** 3 - 3 * (1 - z)), 9, (1, 10), 25, 90, None, lambda v: 1 - v),
    (lambda z: 2 * z / (1 + z ** 2), 2, (1, 5), 30, 110, None, lambda v: v),
    (lambda z: 2 * mpmath.sin(z), 2, (1, 5), 30, 90, None, lambda v: v),
    (lambda z: 2 * z * mpmath.exp(z), 2, (1, 10), 30, 110, None, lambda v: v),
    (lambda z: S2() ** (z + 2) - 2, "log2", (-1, 2), 30, 110, (-0.5, -0.125), lambda v: 2 + v),
]


def coeffs_mp(f, n):
    c = mpmath.taylor(f, 0, n)
    pw, cur = [], [mpmath.mpf(0)] * (n + 1)
    cur = list(c)
    for k in range(1, n + 1):
        pw.append(cur)
        nxt = [mpmath.mpf(0)] * (n + 1)
        for i, a in enumerate(cur):
            if a:
                for j in range(n + 1 - i):
                    nxt[i + j] += a * c[j]
        cur = nxt
    return pw


def lam_of(lam):
    return mpmath.log(2) if lam == "log2" else mpmath.mpf(lam)


def half_machine(p):
    f, lam, (a, b), n, _, br, post = p
    z0 = a / b
    c = np.array([float(v) for v in mpmath.taylor(f, 0, n)])
    pw, cur = [], c.copy()
    for k in range(n):
        pw.append(cur)
        cur = np.convolve(cur, c)[:n + 1]
    L = float(lam_of(lam))
    A = np.array([[(1.0 if k == 0 else 0.0) if j == 0 else (pw[k][j + 1] - L * (j == k)) / L ** (j + 1)
                   for k in range(n)] for j in range(n)])
    e1 = np.zeros(n); e1[0] = 1.0
    cs = np.linalg.solve(A, e1)
    s = lambda x: sum(cs[k] * x ** (k + 1) for k in range(n))
    rhs = np.sqrt(L) * s(z0)
    lo, hi = br if br else (z0, 4 * z0)
    return post(brentq(lambda x: s(x) - rhs, lo, hi, xtol=1e-17, rtol=1e-15))


def half_mp(p):
    f, lam, (a, b), _, n, br, post = p
    z0 = mpmath.mpf(a) / b
    pw = coeffs_mp(f, n)
    L = lam_of(lam)
    A = mpmath.matrix(n, n)
    for j in range(n):
        for k in range(n):
            A[j, k] = (1 if k == 0 else 0) if j == 0 else (pw[k][j + 1] - (L if j == k else 0)) / L ** (j + 1)
    cs = mpmath.lu_solve(A, mpmath.matrix([1] + [0] * (n - 1)))
    s = lambda x: sum(cs[k] * x ** (k + 1) for k in range(n))
    rhs = mpmath.sqrt(L) * s(z0)
    return post(mpmath.findroot(lambda x: s(x) - rhs, z0))


def robbins_machine(n):
    from math import factorial, log
    A = np.array([[(k * log(2)) ** j / factorial(j) - (j == k) for k in range(1, n + 1)] for j in range(n)])
    b = np.zeros(n); b[0] = 1.0
    v = np.linalg.solve(A, b)
    p = np.polynomial.Polynomial(np.concatenate([[-0.5], v]))
    return brentq(p, 0.4, 0.7, xtol=1e-17, rtol=1e-15)


def robbins_mp(n):
    L = mpmath.log(2)
    A = mpmath.matrix([[(k * L) ** j / mpmath.factorial(j) - (j == k) for k in range(1, n + 1)] for j in range(n)])
    v = mpmath.lu_solve(A, mpmath.matrix([1] + [0] * (n - 1)))
    return mpmath.findroot(lambda x: -mpmath.mpf(1) / 2 + sum(v[k] * x ** (k + 1) for k in range(n)), mpmath.mpf(1) / 2)


for k, p in enumerate(PROBS, 1):
    vrun("hd20.%d machine" % k, lambda: half_machine(p), 10, cap=60)
    vrun("hd20.%d wp50" % k, lambda: half_mp(p), 48, cap=120, dps=60)
vrun("hd20.10 machine", lambda: robbins_machine(16), 10, cap=60)
vrun("hd20.10 wp50", lambda: robbins_mp(16), 48, cap=60, dps=60)
