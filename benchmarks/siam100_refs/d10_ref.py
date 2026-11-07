"""Reference for d10, method B (independent of d10_ref.wl's Poisson solve):
power iteration of the Markov operator  (U f)(r) = (f(2-1/r) + f(2+1/r))/2
on [1,3], starting from f = log r.  U^k log r -> lambda (a constant) for every
r, because the ratio chain r_n = x_{n+1}/x_n forgets its start; exp(lambda) is
the a.s. limit of |x_n|^(1/n).  f is held by its values at the N Chebyshev
points of the SECOND kind and pushed through barycentric interpolation
(the WL method uses first-kind points and a direct linear solve).
Usage: python d10_ref.py N dps iters"""
import sys
from mpmath import mp, mpf, cos, pi, log, exp, fdot, nstr
N = int(sys.argv[1]); mp.dps = int(sys.argv[2]); iters = int(sys.argv[3])
xs = [cos(pi * j / (N - 1)) for j in range(N)]          # in [-1,1], r = 2 + x
wts = [(-1) ** j * (mpf(1) / 2 if j in (0, N - 1) else 1) for j in range(N)]
def bary_row(z):
    """Row of the barycentric interpolation matrix at point z."""
    d = [z - x for x in xs]
    for j, dj in enumerate(d):
        if dj == 0:
            return [mpf(1) if k == j else mpf(0) for k in range(N)]
    q = [w / dj for w, dj in zip(wts, d)]
    s = sum(q)
    return [v / s for v in q]
M = []
for x in xs:
    r = 2 + x
    ra, rb = bary_row(-1 / r), bary_row(1 / r)
    M.append([(a + b) / 2 for a, b in zip(ra, rb)])
f = [log(2 + x) for x in xs]
prev = None
for k in range(iters):
    f = [fdot(row, f) for row in M]
    if k % 20 == 19:
        lo, hi = min(f), max(f)
        est = exp((lo + hi) / 2)
        print(k + 1, "spread", nstr(hi - lo, 3), nstr(est, 70), flush=True)
