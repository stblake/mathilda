"""Holdout-A D10 family (experiment 215) references, method A (generalises
d10_ref.py): power iteration of the Markov operator of the ratio chain
r_n = x_{n+1}/x_n,  (U f)(r) = Sum_i p_i f(c_i + s_i/r),  on an interval J
that every map sends into itself, from f = log r.  U^k log r -> lambda, and
exp(lambda) is the a.s. limit of |x_n|^(1/n).  f is held by its values at
N Chebyshev points of the SECOND kind on J (barycentric interpolation).
Method B: hd10_ref.wl (first-kind collocation of the Poisson equation).
Usage: python hd10_ref.py <case> N dps iters"""
import sys
from mpmath import mp, mpf, cos, pi, log, exp, fdot, nstr
mp.dps = int(sys.argv[3])      # BEFORE the probabilities below are formed (1/3 is inexact)
# case: (maps [(c, s, p)], J)
CASES = {
  2: ([(3, -1, mpf(1)/2), (3, 1, mpf(1)/2)], (2, 4)),
  3: ([(4, -1, mpf(1)/2), (4, 1, mpf(1)/2)], (3, 5)),
  4: ([(mpf(5)/2, -1, mpf(1)/2), (mpf(5)/2, 1, mpf(1)/2)], (mpf(3)/2, mpf(7)/2)),
  5: ([(3, -2, mpf(1)/2), (3, 2, mpf(1)/2)], (2, 4)),
  6: ([(1, 1, mpf(1)/2), (2, 1, mpf(1)/2)], (1, 3)),
  7: ([(3, -1, mpf(1)/3), (3, 0, mpf(1)/3), (3, 1, mpf(1)/3)], (2, 4)),
  8: ([(3, -1, mpf(1)/4), (3, 1, mpf(3)/4)], (2, 4)),
  9: ([(1, 1, mpf(1)/3), (2, 1, mpf(1)/3), (3, 1, mpf(1)/3)], (1, 4)),
  10: ([(2, -1, mpf(1)/4), (2, 1, mpf(3)/4)], (1, 3)),
}
case = int(sys.argv[1]); N = int(sys.argv[2]); iters = int(sys.argv[4])
maps, (lo, hi) = CASES[case]
lo, hi = mpf(lo), mpf(hi); mid, half = (lo + hi) / 2, (hi - lo) / 2
xs = [cos(pi * j / (N - 1)) for j in range(N)]
wts = [(-1) ** j * (mpf(1) / 2 if j in (0, N - 1) else 1) for j in range(N)]
def bary_row(r):
    z = (r - mid) / half
    d = [z - x for x in xs]
    for j, dj in enumerate(d):
        if dj == 0:
            return [mpf(1) if k == j else mpf(0) for k in range(N)]
    q = [w / dj for w, dj in zip(wts, d)]
    s = sum(q)
    return [v / s for v in q]
M = []
for x in xs:
    r = mid + half * x
    row = [mpf(0)] * N
    for c, s, p in maps:
        img = c + mpf(s) / r
        assert lo <= img <= hi
        row = [a + p * b for a, b in zip(row, bary_row(img))]
    M.append(row)
f = [log(mid + half * x) for x in xs]
for k in range(iters):
    f = [fdot(row, f) for row in M]
    if k % 20 == 19 or k == iters - 1:
        a, b = min(f), max(f)
        print(case, k + 1, "spread", nstr(b - a, 3), nstr(exp((a + b) / 2), 66), flush=True)
