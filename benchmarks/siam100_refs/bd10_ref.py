# Holdout B, family D10 -- method A: power iteration of the Markov operator
# (U f)(r) = sum_i p_i f(a_i + b_i/r) on f_0 = log r, f_k -> lambda, with f_k
# held by barycentric interpolation at N Chebyshev points of the SECOND kind
# on the invariant interval [lo, hi] (found by iterating the extreme maps),
# mpmath dps 80.  Method B: bd10_ref.gp (Poisson-equation collocation at
# Chebyshev points of the FIRST kind, linear solve).  bd10.1: lambda = mean log c.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bd10_cases import CASES
import mpmath as mp
mp.mp.dps = 80

def interval(cs):
    lo, hi = mp.mpf(1), mp.mpf(1)
    for _ in range(3000):
        ims = [a + mp.mpf(b) / t for (_, a, b) in cs for t in (lo, hi)]
        nlo, nhi = min(ims), max(ims)
        if abs(nlo - lo) + abs(nhi - hi) < mp.mpf(10) ** -75:
            break
        lo, hi = nlo, nhi
    return lo, hi

def lam(cs, N, iters):
    lo, hi = interval(cs)
    if hi - lo < mp.mpf(10) ** -60:
        return mp.log(lo)
    th = [mp.pi * j / (N - 1) for j in range(N)]
    xs = [(lo + hi) / 2 + (hi - lo) / 2 * mp.cos(t) for t in th]
    w = [(-1) ** j * (mp.mpf(1) / 2 if j in (0, N - 1) else 1) for j in range(N)]
    def interp(fv, x):
        num = den = 0
        for j in range(N):
            d = x - xs[j]
            if d == 0:
                return fv[j]
            c = w[j] / d
            num += c * fv[j]; den += c
        return num / den
    pts = [[(mp.mpf(p.numerator) / p.denominator, a + mp.mpf(b) / x) for (p, a, b) in cs] for x in xs]
    f = [mp.log(x) for x in xs]
    prev = None
    for k in range(iters):
        f = [mp.fsum(p * interp(f, y) for (p, y) in row) for row in pts]
        spread = max(f) - min(f)
        if spread < mp.mpf(10) ** -72:
            break
    return (max(f) + min(f)) / 2, spread, k, lo, hi

N = int(sys.argv[1]) if len(sys.argv) > 1 else 90
for lab, cs in CASES.items():
    if all(b == 0 for (_, _, b) in cs):
        L = mp.fsum(mp.mpf(p.numerator) / p.denominator * mp.log(a) for (p, a, _) in cs)
        print(lab, mp.nstr(mp.exp(L), 70), "closed form", mp.nstr(mp.root(24, 4), 70)); continue
    L, spread, k, lo, hi = lam(cs, N, 400)
    print(lab, mp.nstr(mp.exp(L), 70), "N=%d iters=%d spread=%s lo=%s hi=%s" % (N, k, mp.nstr(spread, 3), mp.nstr(lo, 8), mp.nstr(hi, 8)))
    sys.stdout.flush()
