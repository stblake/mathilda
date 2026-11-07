"""Experiment 195 -- SIAM 100-Digit Challenge, App. D Problem 20 (Python column).

No built-in Kneser half-iterate: d20.machine / d20.wp50 SKIP.  The finite
proxy (Robbins' truncated Abel system, see half_exponential.m):
numpy.linalg.solve + scipy.optimize.brentq at machine precision; mpmath
lu_solve + findroot at dps 60 (mpmath has no exact rational solver; the
matrix is built from exact integers/factorials and solved at 60 digits)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
from math import factorial
import numpy as np
import mpmath
from scipy.optimize import brentq

require(["numpy.linalg:solve", "scipy.optimize:brentq", "mpmath:lu_solve", "mpmath:findroot"])

skip("d20.machine", "KneserHalfIterate")
skip("d20.wp50", "KneserHalfIterate")


def h0_machine(n):
    A = np.array([[k ** j / factorial(j) - (j == k) for k in range(1, n + 1)]
                  for j in range(n)])
    b = np.zeros(n); b[0] = 1.0
    v = np.linalg.solve(A, b)
    p = np.polynomial.Polynomial(np.concatenate([[-0.5], v]))   # slog_n(x) + 1/2
    return brentq(p, 0.3, 0.7, xtol=1e-17, rtol=1e-15)


def h0_mp(n):
    A = mpmath.matrix([[mpmath.mpf(k) ** j / mpmath.factorial(j) - (j == k)
                        for k in range(1, n + 1)] for j in range(n)])
    b = mpmath.matrix([1] + [0] * (n - 1))
    v = mpmath.lu_solve(A, b)
    return mpmath.findroot(lambda x: -mpmath.mpf(1) / 2 + mpmath.polyval([v[k] for k in range(n - 1, -1, -1)] + [0], x),
                           mpmath.mpf(1) / 2)


vrun("d20.machine robbins n=20", lambda: h0_machine(20), 10, 3)
vrun("d20.wp50 robbins n=60", lambda: h0_mp(60), 48, dps=60)
