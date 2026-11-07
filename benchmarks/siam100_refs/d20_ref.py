"""Reference A for SIAM100 App. D problem 20 (exp 195): h(h(z)) = exp(z), h(0)?

h(0) for the Kneser half-iterate: h = sexp(slog(z) + 1/2), slog(0) = -1, so
h(0) = sexp(-1/2) = the x with slog(x) = -1/2.

Method A (this file): Robbins' truncated Abel system ("intuitive slog").
  slog_n(x) = -1 + sum_{k=1}^n v_k x^k,  with v solving the n x n system
      sum_{k=1}^n v_k (k^j / j! - [j == k]) = [j == 0],   j = 0..n-1
  (coefficients of x^j in slog_n(e^x) - slog_n(x) = 1).  Solved EXACTLY
  (fractions.Fraction Gaussian elimination), then the root of
  slog_n(x) = -1/2 near 0.4986 by mpmath.findroot.  Robbins' slog_n is
  widely believed (not proved here) to converge to Kneser's slog as n grows;
  the sequence h_n(0) is printed for increasing n to show the digits settle.
  The same routine gives the TRUE values of the finite proxy rows of exp 195
  (fixed n), which are exact-arithmetic answers, not limits.
  h_n(0), n = 80..160: 0.4985633095, 0.4985632904, 0.4985632855, 0.4985632851, 0.4985632851.
Usage: python3.11 d20_ref.py n1 n2 ...
"""
import sys
from fractions import Fraction
from math import factorial
import mpmath as mp

mp.mp.dps = 80


def robbins(n):
    A = [[Fraction(k ** j, factorial(j)) - (1 if j == k else 0) for k in range(1, n + 1)]
         for j in range(n)]
    b = [Fraction(1 if j == 0 else 0) for j in range(n)]
    # Gaussian elimination, exact
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]; b[c], b[piv] = b[piv], b[c]
        inv = 1 / A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] * inv
            if f:
                row_r, row_c = A[r], A[c]
                for k in range(c, n):
                    row_r[k] -= f * row_c[k]
                b[r] -= f * b[c]
    v = [Fraction(0)] * n
    for c in range(n - 1, -1, -1):
        s = b[c] - sum(A[c][k] * v[k] for k in range(c + 1, n))
        v[c] = s / A[c][c]
    return v


def h0(n):
    v = [mp.mpf(q.numerator) / q.denominator for q in robbins(n)]
    f = lambda x: -1 + sum(v[k] * x ** (k + 1) for k in range(n)) + mp.mpf(1) / 2
    return mp.findroot(f, mp.mpf("0.4986"))


if __name__ == "__main__":
    for a in sys.argv[1:]:
        print(a, mp.nstr(h0(int(a)), 70), flush=True)
