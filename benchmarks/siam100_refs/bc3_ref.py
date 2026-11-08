"""Holdout B, family C3 (experiment 228) -- reference values, method B (gmpy2).

Independent of bc3_ref.gp (PARI): entries built here from exact Fractions /
exact phase reduction, arithmetic in gmpy2 mpfr/mpc at 360 bits (~108
digits) in numpy object arrays, and the OTHER Gram matrix where it matters:
  * well-separated spectra: power iteration x <- A (A^H x) on the left
    singular vector (PARI iterates on the right one), to a 1e-100 stable
    Rayleigh quotient;
  * clustered spectra (bc3.3, 3.4, 3.5, 3.7): Rayleigh-quotient iteration on
    G = A A^H (PARI: fixed-shift inverse iteration on A^H A), each step a
    partial-pivot Gaussian elimination written here, started from the shift
    (double estimate)^2.
Prints sigma_max to 100 digits for each row; compare with bc3_ref.gp.
"""
import math
from fractions import Fraction as Fr
import gmpy2
from gmpy2 import mpfr, mpc
import numpy as np

gmpy2.get_context().precision = 360
PI = gmpy2.const_pi()


def ent(k, i, j, n):
    if k == 1: return mpfr(Fr(1, i + j - 1))
    if k == 2: return mpfr(Fr(1, math.factorial(i + j - 1)))
    if k == 3: return mpfr(Fr(1, 1 + (i - j) ** 2))
    if k == 4: return mpfr(Fr(1, i)) if j <= i else mpfr(0)
    if k == 5: return mpfr(Fr(1, 2 ** (i - j))) if i >= j else mpfr(Fr(-1, 3) ** (j - i))
    if k == 6: return mpfr(Fr(1, (i + j - 1) * (i + j - 2) // 2 + j))
    if k == 7: return gmpy2.exp(mpc(0, 2 * PI * ((i * j) % 7) / 7)) / (i + j)
    if k == 8: return mpfr(n + 1 - max(i, j)) if j >= i - 1 else mpfr(0)
    if k == 9: return mpfr(Fr(1, 2 ** i + 3 ** j))
    if k == 10: return 1 / gmpy2.sqrt(mpfr(i + j - 1))


def build(k, m, n):
    return np.array([[ent(k, i, j, n) for j in range(1, n + 1)] for i in range(1, m + 1)],
                    dtype=object)


MPC = type(mpc(0))


def cj(v):
    return v.conjugate() if isinstance(v, MPC) else v


def conjT(A):
    return np.array([[cj(x) for x in row] for row in A.T], dtype=object)


def n2(x):
    return sum((abs(v) ** 2 for v in x), mpfr(0))


def power(A):
    AH = conjT(A)
    x = np.array([mpfr(1)] * A.shape[0], dtype=object)
    l0 = mpfr(0)
    for it in range(100000):
        y = AH.dot(x)
        z = A.dot(y)
        l1 = n2(y) / n2(x)
        x = z / gmpy2.sqrt(n2(z))
        if it > 5 and abs(l1 - l0) < mpfr("1e-100") * l1:
            return gmpy2.sqrt(l1), it
        l0 = l1


def solve(G, b):
    """Gaussian elimination with partial pivoting (object arrays)."""
    n = len(b)
    M = np.concatenate([G.copy(), b.reshape(n, 1)], axis=1)
    for c in range(n):
        p = c + int(np.argmax([abs(v) for v in M[c:, c]]))
        if p != c:
            M[[c, p]] = M[[p, c]]
        piv = M[c, c]
        l = M[c + 1:, c] / piv
        M[c + 1:, c:] = M[c + 1:, c:] - np.outer(l, M[c, c:])
    x = np.empty(n, dtype=object)
    for r in range(n - 1, -1, -1):
        x[r] = (M[r, n] - M[r, r + 1:n].dot(x[r + 1:])) / M[r, r]
    return x


def rqi(A, s0):
    G = A.dot(conjT(A))                      # A A^H
    n = G.shape[0]
    x = np.array([mpfr(1)] * n, dtype=object)
    x = x / gmpy2.sqrt(n2(x))
    mu = mpfr(s0) ** 2
    I = np.array([[mpfr(1) if r == c else mpfr(0) for c in range(n)] for r in range(n)], dtype=object)
    for it in range(30):
        y = solve(G - mu * I, x)
        x = y / gmpy2.sqrt(n2(y))
        xh = np.array([cj(v) for v in x], dtype=object)
        new = xh.dot(G.dot(x))
        new = new.real if isinstance(new, MPC) else new
        if abs(new - mu) < mpfr("1e-100") * new:
            return gmpy2.sqrt(new), it
        mu = new
    raise RuntimeError("RQI did not converge")


JOBS = [(1, "m", 500, 500, None), (1, "w", 60, 60, None), (2, "m", 40, 40, None),
        (3, "m", 400, 400, 3.135326220838264), (3, "w", 60, 60, 3.038012952371806),
        (4, "m", 500, 500, 1.7192159353130747), (4, "w", 80, 80, 1.6144112097670473),
        (5, "m", 400, 400, 1.7499030315692388), (5, "w", 60, 60, 1.7459121257721564),
        (6, "m", 500, 150, None), (6, "w", 100, 30, None),
        (7, "m", 300, 300, 1.07493889956485), (7, "w", 50, 50, 1.007912280679732),
        (8, "m", 300, 300, None), (8, "w", 60, 60, None), (9, "m", 150, 150, None),
        (10, "m", 1000, 1000, None), (10, "w", 120, 120, None)]

if __name__ == "__main__":
    import sys
    only = set(sys.argv[1:])
    for k, tag, m, n, s0 in JOBS:
        if only and ("%d%s" % (k, tag)) not in only:
            continue
        A = build(k, m, n)
        s, it = power(A) if s0 is None else rqi(A, s0)
        print("bc3.%d %s %d %d it=%d %s" % (k, tag, m, n, it, gmpy2.mpfr(s).__format__(".100g")), flush=True)
