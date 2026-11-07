"""Holdout B, family D2 reference, method A (mpmath, dps 90).

Every digit-length series is grouped by the number of base-b digits k of n:
the block of k-digit n contributes  w(k) * H_b(k),
  H_b(k) = sum_{b^(k-1) <= n < b^k} 1/n = psi(b^k) - psi(b^(k-1)) = log b + e_b(k),
with e_b(k) = O(b^-k), so every correction series below converges
geometrically; the log b part is summed in closed form (zeta(2), e - 1,
log 2, or the self-similar linear equation for S).

  bd2.1  sum 1/(n d(n)^2)            = log10 zeta(2) + sum e10(k)/k^2
  bd2.2  sum 1/(n b(n)^2)            = log2 zeta(2) + sum e2(k)/k^2
  bd2.3  sum_{n>=2} 1/(n log^2 n)    : direct to n = N, Euler-Maclaurin tail, exact integral
  bd2.4  sum_{n>=3} 1/(n log n loglog^2 n) : same
  bd2.5  g(1)=1, g(n)=n g(floor(log2 n)): S = (1 + sum_{m>=1} e2(m+1)/g(m))/(1 - log 2)
  bd2.6  f(1)=1, f(n)=n f(d(n))^2:  S = H(9) + log10 T + sum_{k>=2} e10(k)/f(k)^2,
         T = sum_{k>=2} 1/f(k)^2 = sum_{k=2}^9 1/k^2 + sum_{m>=2} (psi1(10^(m-1)) - psi1(10^m))/f(m)^4
  bd2.7  sum 1/(n d(n) d(d(n))^2)    = log10 * [bd2.1] + sum_k e10(k)/(k d(k)^2)
  bd2.8  sum 1/(n d(n)!)             = log10 (e - 1) + sum e10(k)/k!
  bd2.9  sum (-1)^b(n)/(n b(n))      = -log2^2 + sum (-1)^k e2(k)/k
  bd2.10 f(1)=1, f(2)=2, f(n)=2n f(t(n)) (t = ternary length):
         S = (3/2 - log3/2 + sum_{k>=2} e3(k)/(2 f(k)))/(1 - log3/2)
Method B (independent code and decomposition): bd2_ref.gp."""
import mpmath as mp
mp.mp.dps = 90
KMAX = 330


def H(b, k):
    return mp.digamma(mp.mpf(b) ** k) - mp.digamma(mp.mpf(b) ** (k - 1))


def e(b, k):
    return H(b, k) - mp.log(b)


def ndig(n, b=10):
    c = 0
    while n:
        n //= b
        c += 1
    return c


out = {}
out["bd2.1"] = mp.log(10) * mp.zeta(2) + mp.fsum(e(10, k) / k ** 2 for k in range(1, KMAX))
out["bd2.2"] = mp.log(2) * mp.zeta(2) + mp.fsum(e(2, k) / k ** 2 for k in range(1, KMAX))
# smooth series: direct head to N, then Euler-Maclaurin with the exact tail
# integral (1/log N, 1/loglog N) and derivatives at N from mpmath.taylor
# (numerical differentiation at raised precision; method B uses exact
# power-series derivatives instead).  (mpmath.sumem's own quadrature of the
# integral to infinity is NOT used: on this logarithmic tail it is wrong in
# the 3rd digit.)
N = 2000


def em(f, a, I):
    S = mp.fsum(f(mp.mpf(n)) for n in range(a, N)) + I + f(mp.mpf(N)) / 2
    c = mp.taylor(f, mp.mpf(N), 61)
    for k in range(1, 31):
        S -= mp.bernoulli(2 * k) / (2 * k) * c[2 * k - 1]
    return S


f3 = lambda x: 1 / (x * mp.log(x) ** 2)
out["bd2.3"] = em(f3, 2, 1 / mp.log(N))
f4 = lambda x: 1 / (x * mp.log(x) * mp.log(mp.log(x)) ** 2)
out["bd2.4"] = em(f4, 3, 1 / mp.log(mp.log(N)))


def g5(n):
    return 1 if n == 1 else n * g5(n.bit_length() - 1)


E = mp.fsum(e(2, m + 1) / g5(m) for m in range(1, KMAX))
out["bd2.5"] = (1 + E) / (1 - mp.log(2))


def f6(n):
    return 1 if n == 1 else n * f6(ndig(n)) ** 2


T = mp.fsum(mp.mpf(1) / k ** 2 for k in range(2, 10)) + mp.fsum(
    (mp.psi(1, mp.mpf(10) ** (m - 1)) - mp.psi(1, mp.mpf(10) ** m)) / mp.mpf(f6(m)) ** 4
    for m in range(2, 120))
out["bd2.6"] = (mp.fsum(mp.mpf(1) / n for n in range(1, 10)) + mp.log(10) * T
                + mp.fsum(e(10, k) / mp.mpf(f6(k)) ** 2 for k in range(2, KMAX)))
out["bd2.7"] = mp.log(10) * out["bd2.1"] + mp.fsum(
    e(10, k) / (k * ndig(k) ** 2) for k in range(1, KMAX))
out["bd2.8"] = mp.log(10) * (mp.e - 1) + mp.fsum(e(10, k) / mp.factorial(k) for k in range(1, KMAX))
out["bd2.9"] = -mp.log(2) ** 2 + mp.fsum((-1) ** k * e(2, k) / k for k in range(1, KMAX))


def f10(n):
    return n if n <= 2 else 2 * n * f10(ndig(n, 3))


E = mp.fsum(e(3, k) / (2 * f10(k)) for k in range(2, KMAX))
out["bd2.10"] = (mp.mpf(3) / 2 - mp.log(3) / 2 + E) / (1 - mp.log(3) / 2)
for k, v in out.items():
    print(k, mp.nstr(v, 78))
