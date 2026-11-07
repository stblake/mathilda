"""Experiment 212 -- SIAM 100-digit HOLDOUT A, family D7 (Python column).

The 1-D statements of region_areas.m: ends by root-finding from the same
starting values (scipy.optimize.newton, secant mode; mpmath.findroot), then
the chord length integrated with the ends and kinks as breakpoints
(scipy.integrate.quad, epsabs 1e-14, epsrel 1e-13; mpmath.quad at 60 digits).
numpy.cbrt is the real cube root; the real odd root in mpmath (and the fifth
root in numpy) is written sign(t) |t|^(1/n)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import numpy as np
import mpmath as M
from scipy.optimize import newton
from scipy.integrate import quad

require(["scipy.optimize:newton", "scipy.integrate:quad", "mpmath:findroot", "mpmath:quad"])

Q = lambda f, a, b: quad(f, a, b, epsabs=1e-14, epsrel=1e-13, limit=200)[0]
# machine spellings
g16 = lambda y: y * y / 16 + np.exp(np.sin(y))
g32 = lambda y: y * y / 32 + np.exp(np.sin(y))
g100 = lambda y: y * y / 100 + np.exp(np.sin(y))
root5 = lambda t: np.sign(t) * abs(t) ** 0.2
def fg_m(g, inv, sc, islands):
    tot = 0.0
    for lo0, hi0, ks0 in islands:
        a = 0.0 if lo0 == 0 else newton(lambda y: g(y) - 1, lo0, tol=1e-15)
        b = 0.0 if hi0 == 0 else newton(lambda y: g(y) - 1, hi0, tol=1e-15)
        ks = [newton(lambda y: -sc * np.log(g(y)) - y * y, k0, tol=1e-15) for k0 in ks0]
        def w(y):
            s = np.sqrt(-sc * np.log(g(y)))
            return inv(-y + s) - inv(-y - s)
        tot += quad(w, a, b, points=ks or None, epsabs=1e-14, epsrel=1e-13, limit=200)[0]
    return tot
# 60-digit spellings
G16 = lambda y: y ** 2 / 16 + M.exp(M.sin(y))
G32 = lambda y: y ** 2 / 32 + M.exp(M.sin(y))
G100 = lambda y: y ** 2 / 100 + M.exp(M.sin(y))
rc3 = lambda t: M.sign(t) * M.cbrt(abs(t))
rc5 = lambda t: M.sign(t) * M.root(abs(t), 5)
def fg_50(g, inv, sc, islands):
    tot = 0
    for lo0, hi0, ks0 in islands:
        a = 0 if lo0 == 0 else M.findroot(lambda y: g(y) - 1, M.mpf(lo0))
        b = 0 if hi0 == 0 else M.findroot(lambda y: g(y) - 1, M.mpf(hi0))
        ks = [M.findroot(lambda y: -sc * M.log(g(y)) - y ** 2, M.mpf(k0)) for k0 in ks0]
        def w(y):
            s = M.sqrt(-sc * M.log(g(y)))
            return inv(-y + s) - inv(-y - s)
        tot += M.quad(w, [a] + ks + [b])
    return tot
I3 = [("-2.6", 0, ["-0.8"])]
I4 = [("-2.8", 0, []), ("3.8", "4.45", [])]
I5 = [("-2.8", 0, ["-0.8"]), ("3.8", "4.45", [])]
I6 = [("-2.8", 0, ["-1.3"]), ("3.8", "4.45", [])]
I10 = [("-7.94", "-7.04", []), ("-3.04", 0, ["-0.86"]), ("3.25", "5.85", [])]
fl = lambda I: [(0 if a == 0 else float(a), 0 if b == 0 else float(b), [float(k) for k in ks]) for a, b, ks in I]
L8 = math.log(2 + math.sqrt(3))
def a2_m():
    a = newton(lambda x: math.cos(x) - x * x, 0.8, tol=1e-15)
    return Q(lambda x: math.cos(x) - x * x, -a, a)
def a2_50():
    a = M.findroot(lambda x: M.cos(x) - x ** 2, M.mpf("0.8"))
    return M.quad(lambda x: M.cos(x) - x ** 2, [-a, a])
def a9_m():
    a = newton(lambda x: 3 - x * x - math.exp(x), -1.7, tol=1e-15)
    b = newton(lambda x: 3 - x * x - math.exp(x), 0.8, tol=1e-15)
    return Q(lambda x: 3 - x * x - math.exp(x), a, b)
def a9_50():
    a = M.findroot(lambda x: 3 - x ** 2 - M.exp(x), M.mpf("-1.7"))
    b = M.findroot(lambda x: 3 - x ** 2 - M.exp(x), M.mpf("0.8"))
    return M.quad(lambda x: 3 - x ** 2 - M.exp(x), [a, b])
CASES = {
  1: (lambda: Q(lambda x: 2 * (1 - x ** 4) ** 0.25, -1, 1),
      lambda: M.quad(lambda x: 2 * (1 - x ** 4) ** (M.mpf(1) / 4), [-1, 1])),
  2: (a2_m, a2_50),
  3: (lambda: fg_m(g16, np.cbrt, 1, fl(I3)), lambda: fg_50(G16, rc3, 1, I3)),
  4: (lambda: fg_m(g32, np.arcsinh, 1, fl(I4)), lambda: fg_50(G32, M.asinh, 1, I4)),
  5: (lambda: fg_m(g32, root5, 1, fl(I5)), lambda: fg_50(G32, rc5, 1, I5)),
  6: (lambda: fg_m(g32, np.cbrt, 2, fl(I6)), lambda: fg_50(G32, rc3, 2, I6)),
  7: (lambda: Q(lambda x: math.pi - 2 * math.asin(1 - math.sin(x)), 0, math.pi),
      lambda: M.quad(lambda x: M.pi - 2 * M.asin(1 - M.sin(x)), [0, M.pi])),
  8: (lambda: Q(lambda x: math.log(4 - math.exp(x)) + x, -L8, L8),
      lambda: M.quad(lambda x: M.log(4 - M.exp(x)) + x, [-M.log(2 + M.sqrt(3)), M.log(2 + M.sqrt(3))])),
  9: (a9_m, a9_50),
  10: (lambda: fg_m(g100, np.cbrt, 1, fl(I10)), lambda: fg_50(G100, rc3, 1, I10)),
}
for k in range(1, 11):
    fm, f5 = CASES[k]
    vrun("hd7.%d machine" % k, fm, 10, 3, cap=30)
    vrun("hd7.%d wp50" % k, f5, 48, cap=60, dps=60)
