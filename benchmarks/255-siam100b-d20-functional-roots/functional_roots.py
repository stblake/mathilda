"""Experiment 255 -- SIAM 100-digit SEALED HOLDOUT B, family D20 (Python column).

See functional_roots.m for the ten problems.  Machine rows: numpy.linalg.solve
+ scipy.optimize.brentq (Robbins rows), scipy.linalg.sqrtm (Carleman rows).
wp50 rows: mpmath lu_solve + findroot, mpmath.sqrtm, all at dps 60 (mpmath
has no exact rational solver; the matrices are built from exact integers /
factorials / log(b) and solved at 60 digits)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
from math import factorial, comb
import numpy as np
import scipy.linalg as sl
from scipy.optimize import brentq
from scipy.special import stirling2
import mpmath

require(["numpy.linalg:solve", "scipy.optimize:brentq", "scipy.linalg:sqrtm", "scipy.special:stirling2",
         "mpmath:lu_solve", "mpmath:findroot", "mpmath:sqrtm"])

LOGS = {"e": (1.0, lambda: mpmath.mpf(1)), "2": (np.log(2.0), lambda: mpmath.log(2)),
        "sqrt2": (np.log(2.0) / 2, lambda: mpmath.log(2) / 2), "10": (np.log(10.0), lambda: mpmath.log(10))}


def robbins_np(n, base):
    L = LOGS[base][0]
    A = np.array([[(k * L) ** j / factorial(j) - (j == k) for k in range(1, n + 1)] for j in range(n)])
    b = np.zeros(n); b[0] = 1.0
    return np.linalg.solve(A, b)


def robbins_mp(n, base):
    L = LOGS[base][1]()
    A = mpmath.matrix([[(k * L) ** j / mpmath.factorial(j) - (j == k) for k in range(1, n + 1)]
                       for j in range(n)])
    v = mpmath.lu_solve(A, mpmath.matrix([1] + [0] * (n - 1)))
    return [v[i] for i in range(n)]


def root_machine(n, base, level, a, b):
    v = robbins_np(n, base)
    p = np.polynomial.Polynomial(np.concatenate([[-1.0 - level], v]))   # slog_n(x) - level
    return brentq(p, a, b, xtol=1e-17, rtol=1e-15)


def root_mp(n, base, level, x0):
    v = robbins_mp(n, base)
    return mpmath.findroot(lambda x: -1 - mpmath.mpf(level) + mpmath.polyval(v[::-1] + [0], x), mpmath.mpf(x0))


c7 = lambda j, k: comb(j, k - j) * 2 ** (2 * j - k) if j <= k <= 2 * j else 0
c8 = lambda j, k: factorial(j) * int(stirling2(k, j, exact=True)) / mpmath.mpf(factorial(k)) if k >= j else 0
c9 = lambda j, k: 2 ** j * comb(j + (k - j) // 2 - 1, (k - j) // 2) if (k >= j and (k - j) % 2 == 0) else 0
c8f = lambda j, k: factorial(j) * float(stirling2(k, j, exact=True)) / factorial(k) if k >= j else 0.0


def half_machine(coef, n, x0):
    C = np.array([[float(coef(j, k)) for k in range(1, n + 1)] for j in range(1, n + 1)])
    S = np.real(sl.sqrtm(C))
    return float(S[0] @ (x0 ** np.arange(1, n + 1)))


def half_mp(coef, n, x0):
    C = mpmath.matrix([[coef(j, k) for k in range(1, n + 1)] for j in range(1, n + 1)])
    S = mpmath.sqrtm(C)
    # the principal root of a real matrix is real: mpmath carries ~1e-60 imaginary noise
    return mpmath.re(sum(S[0, k - 1] * mpmath.mpf(x0) ** k for k in range(1, n + 1)))


R = mpmath.mpf
CASES = [
    ("bd20.1", lambda: root_machine(10, "e", -0.5, 0.4, 0.6), lambda: root_mp(10, "e", R(-1) / 2, "0.5"), 10),
    ("bd20.2", lambda: root_machine(20, "2", -0.5, 0.45, 0.65), lambda: root_mp(20, "2", R(-1) / 2, "0.54"), 10),
    ("bd20.3", lambda: root_machine(30, "sqrt2", -0.5, 0.55, 0.7), lambda: root_mp(30, "sqrt2", R(-1) / 2, "0.63"), 10),
    ("bd20.4", lambda: root_machine(15, "10", -0.5, 0.3, 0.5), lambda: root_mp(15, "10", R(-1) / 2, "0.39"), 10),
    ("bd20.5", lambda: root_machine(30, "e", -0.75, 0.2, 0.35), lambda: root_mp(30, "e", R(-3) / 4, "0.27"), 10),
    ("bd20.6", lambda: robbins_np(24, "e")[0], lambda: robbins_mp(24, "e")[0], 10),
    ("bd20.7", lambda: half_machine(c7, 20, 0.5), lambda: half_mp(c7, 20, R(1) / 2), 10),
    ("bd20.8", lambda: half_machine(c8f, 15, 0.5), lambda: half_mp(c8, 15, R(1) / 2), 10),
    ("bd20.9", lambda: half_machine(c9, 25, 0.3), lambda: half_mp(c9, 25, R(3) / 10), 10),
    ("bd20.10", lambda: root_machine(45, "e", -0.5, 0.4, 0.6), lambda: root_mp(45, "e", R(-1) / 2, "0.5"), 7),
]
for lab, fm, fw, dm in CASES:
    vrun(lab + " machine", fm, dm, 3, cap=60)
    vrun(lab + " wp50", fw, 48, cap=90, dps=60)
