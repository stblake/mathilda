"""Holdout A, family hd2 (experiment 207): self-similar digit-length series.

Each problem: f(n) = seed[n] for n <= c, else f(n) = n^w * L^e * f(L)^q with
L = (number of base-b digits of n) - shift; sum a weighted series of f(n)^-p.

Exact block decomposition: for n of base-b length k (n > c), f(n)^-p =
n^-(w p) * L^-(e p) * f(L)^-(q p), and the block sums of n^-sigma over a
residue class are Hurwitz-zeta / digamma differences, so

  sum_{n >= N} n^-s f(n)^-p  =  direct part  +  sum_k  L^-(ep) f(L)^-(qp) Z_k(sigma).

For sigma > 1, Z_k decays geometrically and the k-sum is summed to convergence.
For sigma == 1, Z_k = log b + delta_k exactly (delta_k = O(b^-k)), so the k-sum
splits into a geometric delta part plus log b times a series of the SAME kind
one level down (exponents (qp, ep)).  When that is the same series, starting
the top level at N with len(b^(N-1)) == N makes the equation exact:
X = D + F + log(b) X.  Usage: hd2_ref.py dps"""
import sys
from functools import lru_cache
from mpmath import mp, mpf, zeta, digamma, log, sqrt, nstr

# (b, c, seed, w, e, q, shift, p, weight)  weight: 'all', 'alt', 'odd'
PROBLEMS = {
    0: (2, 2, {1: 1, 2: 2}, 1, 0, 1, 0, 1, 'all'),   # the original (check)
    1: (2, 2, {1: 1, 2: 2}, 2, 0, 1, 0, 1, 'all'),
    2: (2, 2, {1: 1, 2: 3}, 1, 0, 1, 0, 1, 'all'),
    3: (2, 2, {1: 1, 2: 2}, 1, 0, 1, 0, 1, 'alt'),
    4: (10, 9, None, 1, 1, 1, 0, 1, 'all'),
    5: (3, 2, {1: 1, 2: 2}, 1, 0, 2, 0, 1, 'all'),
    6: (2, 2, {1: 1, 2: 2}, 1, 0, 1, 0, 2, 'all'),
    7: (2, 1, {1: 1}, 1, 0, 1, 1, 1, 'all'),
    8: (10, 9, None, 1, 1, 1, 0, 1, 'alt'),
    9: (2, 2, {1: 1, 2: 2}, 1, mpf(1) / 2, 1, 0, 1, 'all'),
    10: (2, 2, {1: 1, 2: 2}, 1, 0, 1, 0, 1, 'odd'),
}


def setup(P):
    b, c, seed, w, e, q, shift, p, weight = P
    if seed is None:
        seed = {n: n for n in range(1, c + 1)}

    def blen(n):
        k = 0
        while n:
            n //= b
            k += 1
        return k

    @lru_cache(maxsize=None)
    def f(n):
        if n <= c:
            return mpf(seed[n])
        L = blen(n) - shift
        return mpf(n) ** w * mpf(L) ** e * f(L) ** q

    def Lof(k):
        return k - shift

    return b, c, f, blen, Lof, w, e, q


def zsum(sig, a, z, cls):
    """sum of n^-sig for a <= n <= z in class cls ('all', 'odd', 'even')."""
    if cls == 'all':
        if sig == 1:
            return digamma(z + 1) - digamma(a)
        return zeta(sig, a) - zeta(sig, z + 1)
    # n = 2m + r
    r = 1 if cls == 'odd' else 0
    m0 = (a - r + 1) // 2       # smallest m with 2m+r >= a
    m1 = (z - r) // 2
    if sig == 1:
        return (digamma(m1 + 1 + mpf(r) / 2) - digamma(m0 + mpf(r) / 2)) / 2
    return (zeta(sig, m0 + mpf(r) / 2) - zeta(sig, m1 + 1 + mpf(r) / 2)) / mpf(2) ** sig


def solve(P):
    b, c, f, blen, Lof, w, e, q = setup(P)
    p0, weight = P[7], P[8]
    eps = mpf(10) ** (-mp.dps - 5)

    def blocks(p, s, k_from, cls, sign):
        """sum over k >= k_from of [block k sums of n^-s f(n)^-p], tail-free part
        for sigma > 1; for sigma == 1 returns (F, recursion-coefficient)."""
        sig = s + w * p
        tot = mpf(0)
        k = k_from
        while True:
            lo, hi = b ** (k - 1), b ** k - 1
            L = Lof(k)
            g = mpf(L) ** (-e * p) * f(L) ** (-q * p) if L >= 1 else None
            if cls == 'alt':
                Z = zsum(sig, lo, hi, 'odd') - zsum(sig, lo, hi, 'even')
            elif sig == 1:
                Z = zsum(1, lo, hi, cls) - (log(b) if cls == 'all' else log(b) / 2)
            else:
                Z = zsum(sig, lo, hi, cls)
            term = g * Z
            tot += term
            if k > k_from + 5 and abs(term) < eps * (1 + abs(tot)):
                return tot
            k += 1

    @lru_cache(maxsize=None)
    def S(p, s, N):
        """sum_{n >= N} n^-s f(n)^-p over all n, N > c."""
        sig = s + w * p
        k_from = blen(N) + 1          # full blocks from here
        D = sum(mpf(n) ** (-s) * f(n) ** (-p) for n in range(N, b ** (k_from - 1)))
        F = blocks(p, s, k_from, 'all', 1)
        if sig != 1:
            return D + F
        # remaining: log b * sum_{k >= k_from} L^-(ep) f(L)^-(qp)
        p2, s2 = q * p, e * p
        if P[6]:
            # L = k - shift runs over all integers >= k_from - shift
            N2 = k_from - P[6]
        else:
            N2 = k_from
        if (p2, s2) == (p, s) and N2 == N:
            return (D + F) / (1 - log(b))
        return D + F + log(b) * Sfrom(p2, s2, N2)

    def Sfrom(p, s, N):
        # sum_{n >= N} n^-s f(n)^-p allowing N <= c (seed part direct)
        if N > c:
            return S(p, s, N)
        head = sum(mpf(n) ** (-s) * f(n) ** (-p) for n in range(N, c + 1))
        return head + Sfrom(p, s, c + 1)

    # choose a self-consistent top start for the fixed point: N with blen(N)+1 == N (+shift)
    if weight == 'all':
        N = c + 1
        if w * p0 == 1 and (q * p0, e * p0) == (p0, 0):
            while blen(N) + 1 - P[6] != N:
                N += 1
                if N > 64:
                    raise RuntimeError
        head = sum(f(n) ** (-p0) for n in range(1, N))
        return head + S(p0, 0, N)
    # weighted top level (alt / odd): blocks from k0, then recursion into the plain series
    k0 = blen(c + 1) + 3
    lim = b ** (k0 - 1)
    sgn = (lambda n: (-1) ** (n + 1)) if weight == 'alt' else (lambda n: n % 2)
    head = sum(sgn(n) * f(n) ** (-p0) for n in range(1, lim))
    F = blocks(p0, 0, k0, weight, 1)
    if weight == 'alt':
        return head + F
    # odd, sigma = 1: Z_k = log(b)/2 + delta_k ; recursion into plain series over L
    sig = w * p0
    assert sig == 1
    return head + F + log(b) / 2 * Sfrom(q * p0, e * p0, k0 - P[6])


if __name__ == "__main__":
    mp.dps = int(sys.argv[1])
    which = [int(a) for a in sys.argv[2:]] or sorted(PROBLEMS)
    for k in which:
        print(k, nstr(solve(PROBLEMS[k]), mp.dps - 10), flush=True)
