"""Reference A for holdout-B family D20 (exp 255): finite linear-algebra
proxies of functional square roots / Abel functions.

Robbins' truncated Abel system for base b (L = log b):
    slog_n(x) = -1 + sum_{k=1}^n v_k x^k,
    sum_{k=1}^n v_k ((k L)^j / j! - [j == k]) = [j == 0],  j = 0..n-1
(the coefficients of x^j in slog_n(b^x) - slog_n(x) = 1).
  base e: solved EXACTLY (fractions.Fraction Gaussian elimination);
  other bases: mpmath lu_solve at 200 digits;
  then the root by mpmath.findroot at 200 digits.
Carleman square roots (C_jk = [x^k] g(x)^j, j,k = 1..n, upper triangular; the
first row of the principal C^(1/2) is the degree-n Taylor polynomial of the
half-iterate h, h(h(x)) = g(x)), evaluated at x0:
  bd20.7  g = 2x + x^2:     h = (1+x)^sqrt2 - 1, closed form sum binom(sqrt2,k) x0^k
  bd20.8  g = exp(x) - 1:   exact rational Taylor coefficients of the formal
                            half-iterate from h(h(x)) = g(x) solved term by term
  bd20.9  g = 2x/(1-x^2):   h = tan(sqrt2 arctan x); coefficients from the
                            ODE (1+x^2) h' = sqrt2 (1 + h^2), 200 digits
None of method A uses a matrix square root; method B (bd20_ref.gp, and
mpmath.sqrtm for bd20.8) does.
Also prints, per machine row, the digits numpy (LAPACK LU / scipy sqrtm) keeps,
and the input-perturbation sensitivity (entries perturbed by 2^-53 relative).
Usage: python3.11 bd20_ref.py
"""
from fractions import Fraction
from math import factorial, comb
import random
import numpy as np
import scipy.linalg as sl
import mpmath as mp

mp.mp.dps = 200
random.seed(1)


def robbins_exact(n):
    A = [[Fraction(k ** j, factorial(j)) - (1 if j == k else 0) for k in range(1, n + 1)] for j in range(n)]
    b = [Fraction(1 if j == 0 else 0) for j in range(n)]
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]; b[c], b[piv] = b[piv], b[c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            if f:
                for k in range(c, n):
                    A[r][k] -= f * A[c][k]
                b[r] -= f * b[c]
    v = [Fraction(0)] * n
    for c in range(n - 1, -1, -1):
        v[c] = (b[c] - sum(A[c][k] * v[k] for k in range(c + 1, n))) / A[c][c]
    return [mp.mpf(q.numerator) / q.denominator for q in v]


def robbins_matrix(n, L, pert=0):
    A = mp.matrix(n, n)
    for j in range(n):
        for k in range(1, n + 1):
            a = (k * L) ** j / mp.factorial(j) - (1 if j == k else 0)
            if pert:
                a *= 1 + mp.mpf(random.uniform(-1, 1)) * mp.mpf(2) ** -53
            A[j, k - 1] = a
    return A


def robbins_mp(n, L, pert=0):
    A = robbins_matrix(n, L, pert)
    v = mp.lu_solve(A, mp.matrix([1] + [0] * (n - 1)))
    return [v[i] for i in range(n)]


def slog_root(v, level, x0):
    f = lambda x: -1 + sum(v[k] * x ** (k + 1) for k in range(len(v))) - level
    return mp.findroot(f, mp.mpf(x0))


def robbins_np(n, L):
    A = np.array([[float((k * L) ** j / mp.factorial(j)) - (j == k) for k in range(1, n + 1)] for j in range(n)])
    b = np.zeros(n); b[0] = 1
    return list(np.linalg.solve(A, b))


def digits(a, b):
    return float(-mp.log10(abs((a - b) / b))) if a != b else 99


# ---- Robbins rows: (label, base log, n, quantity) ----
ROWS = [("bd20.1", 1, 10, ("root", -mp.mpf(1) / 2, "0.5")),
        ("bd20.2", "2", 20, ("root", -mp.mpf(1) / 2, "0.54")),
        ("bd20.3", "sqrt2", 30, ("root", -mp.mpf(1) / 2, "0.63")),
        ("bd20.4", "10", 15, ("root", -mp.mpf(1) / 2, "0.39")),
        ("bd20.5", 1, 30, ("root", -mp.mpf(3) / 4, "0.27")),
        ("bd20.6", 1, 24, ("v1",)),
        ("bd20.10", 1, 45, ("root", -mp.mpf(1) / 2, "0.4986"))]


def Lof(tag):
    return {1: mp.mpf(1), "2": mp.log(2), "sqrt2": mp.log(2) / 2, "10": mp.log(10)}[tag]


def quantity(v, q):
    if q[0] == "v1":
        return v[0]
    return slog_root(v, q[1], q[2])


for lab, tag, n, q in ROWS:
    L = Lof(tag)
    v = robbins_exact(n) if tag == 1 else robbins_mp(n, L)
    val = quantity(v, q)
    print(lab, "n=%d" % n, mp.nstr(val, 70), flush=True)
    with mp.workdps(30):
        pass
    vnp = robbins_np(n, L)
    vn = quantity([mp.mpf(t) for t in vnp], q)
    sens = min(digits(quantity(robbins_mp(n, L, pert=1), q), val) for _ in range(3))
    print("   numpy digits %.1f, perturbation-sensitivity digits %.1f" % (digits(vn, val), sens), flush=True)


# ---- Carleman rows ----
def carleman_np(coef, n):
    return np.array([[float(coef(j, k)) for k in range(1, n + 1)] for j in range(1, n + 1)])


c7 = lambda j, k: comb(j, k - j) * 2 ** (2 * j - k) if j <= k <= 2 * j else 0
c9 = lambda j, k: (2 ** j * comb(j + (k - j) // 2 - 1, (k - j) // 2)) if (k >= j and (k - j) % 2 == 0) else 0


def stirling2(nn, kk):
    return sum((-1) ** i * comb(kk, i) * (kk - i) ** nn for i in range(kk + 1)) // factorial(kk)


c8 = lambda j, k: Fraction(factorial(j) * stirling2(k, j), factorial(k)) if k >= j else 0


def sqrt_row_np(coef, n, x0):
    S = sl.sqrtm(carleman_np(coef, n))
    return float(np.real(S[0]) @ (x0 ** np.arange(1, n + 1)))


# bd20.7
n, x0 = 20, mp.mpf(1) / 2
s2 = mp.sqrt(2)
v7 = sum(mp.binomial(s2, k) * x0 ** k for k in range(1, n + 1))
print("bd20.7 n=20", mp.nstr(v7, 70), " scipy sqrtm digits %.1f" % digits(sqrt_row_np(c7, n, 0.5), v7), flush=True)

# bd20.8: formal half-iterate of exp(x)-1, exact rationals


def compose(a, b, n):            # a(b(x)) for series with a[0]=b[0]=0, lists index = degree
    res = [Fraction(0)] * (n + 1)
    p = [Fraction(1)] + [Fraction(0)] * n        # b^0
    for m in range(1, n + 1):
        q = [Fraction(0)] * (n + 1)
        for i, pi in enumerate(p):
            if pi:
                for j in range(1, n + 1 - i):
                    if b[j]:
                        q[i + j] += pi * b[j]
        p = q
        if a[m]:
            for d in range(n + 1):
                res[d] += a[m] * p[d]
    return res


n8, x8 = 15, Fraction(1, 2)
g = [Fraction(0)] + [Fraction(1, factorial(k)) for k in range(1, n8 + 1)]
h = [Fraction(0), Fraction(1)] + [Fraction(0)] * (n8 - 1)
for k in range(2, n8 + 1):
    hh = compose(h, h, n8)
    h[k] = (g[k] - hh[k]) / 2        # [x^k] h(h) = 2 h_k + (lower terms)
assert compose(h, h, n8) == g
v8q = sum(h[k] * x8 ** k for k in range(1, n8 + 1))
v8 = mp.mpf(v8q.numerator) / v8q.denominator
print("bd20.8 n=15", mp.nstr(v8, 70), " coefficients h_2..h_6:", [str(h[k]) for k in range(2, 7)], flush=True)
print("   scipy sqrtm digits %.1f" % digits(sqrt_row_np(c8, n8, 0.5), v8), flush=True)
# method B for bd20.8: mpmath sqrtm (Denman-Beavers) of the exact Carleman matrix
with mp.workdps(120):
    C = mp.matrix([[mp.mpf(c8(j, k).numerator) / c8(j, k).denominator if c8(j, k) else 0
                    for k in range(1, n8 + 1)] for j in range(1, n8 + 1)])
    S = mp.sqrtm(C)
    v8b = sum(S[0, k - 1] * mp.mpf(1) / 2 ** k for k in range(1, n8 + 1))
    print("   B mpmath.sqrtm:", mp.nstr(v8b, 70), " agree digits %.1f" % digits(v8b, v8), flush=True)

# bd20.9: tan(sqrt2 arctan x), coefficients from (1+x^2) h' = sqrt2 (1 + h^2)
n9, x9 = 25, mp.mpf(3) / 10
hc = [mp.mpf(0)] * (n9 + 2)
for m in range(0, n9):
    # [x^m] of (1+x^2) h' = (m+1) h_{m+1} + (m-1) h_{m-1};  rhs [x^m] sqrt2 (1 + h^2)
    rhs = s2 * ((1 if m == 0 else 0) + sum(hc[i] * hc[m - i] for i in range(0, m + 1)))
    lhs_rest = (m - 1) * hc[m - 1] if m >= 1 else 0
    hc[m + 1] = (rhs - lhs_rest) / (m + 1)
v9 = sum(hc[k] * x9 ** k for k in range(1, n9 + 1))
print("bd20.9 n=25", mp.nstr(v9, 70), " (full h(0.3) =", mp.nstr(mp.tan(s2 * mp.atan(x9)), 20), ")",
      " scipy sqrtm digits %.1f" % digits(sqrt_row_np(c9, n9, 0.3), v9), flush=True)
