"""Experiment 242 -- holdout B for SIAM Appendix D Problem 7 (Python column).

The 1-D statements of region_areas.m: roots from the same starting values
(scipy.optimize.newton, secant mode; mpmath.findroot at 60 digits), then
the slice length integrated over the same pieces (scipy.integrate.quad, the
interior kink/cusp as a breakpoint; mpmath.quad on the same subintervals).
Two wp50 rows (bd7.3, bd7.9) take the real part of an mpmath result whose
imaginary part is a ~1e-100 rounding artefact (asin/acos at an argument
rounded past 1).  numpy.cbrt is the real cube root; mpmath.cbrt is the principal complex one,
so the wp50 rows take the real cube root as sign(t) * cbrt(|t|)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import numpy as np
import mpmath as M
from scipy.optimize import newton
from scipy.integrate import quad

require(["scipy.optimize:newton", "scipy.integrate:quad", "mpmath:findroot", "mpmath:quad"])

def nt(f, s):
    return newton(f, s, tol=1e-15)
def q(f, a, b, pts=None):
    return quad(f, a, b, points=pts, epsabs=1e-14, epsrel=1e-13, limit=200)[0]
mf = lambda s: M.mpf(s)
rc = lambda t: M.sign(t) * M.cbrt(abs(t))

# 1
def a1m():
    x0 = nt(lambda x: np.cos(x) - x * x, 0.8)
    return q(lambda x: np.cos(x) - x * x, -x0, x0)
def a1w():
    x0 = M.findroot(lambda x: M.cos(x) - x ** 2, mf("0.8"))
    return M.quad(lambda x: M.cos(x) - x ** 2, [-x0, x0])
# 2
def a2m():
    a = nt(lambda x: np.exp(x) + x * x - 3, -1.7); b = nt(lambda x: np.exp(x) + x * x - 3, 0.8)
    return q(lambda x: 3 - x * x - np.exp(x), a, b)
def a2w():
    a = M.findroot(lambda x: M.exp(x) + x ** 2 - 3, mf("-1.7")); b = M.findroot(lambda x: M.exp(x) + x ** 2 - 3, mf("0.8"))
    return M.quad(lambda x: 3 - x ** 2 - M.exp(x), [a, b])
# 3
def a3m():
    return q(lambda x: np.pi - 2 * np.arcsin(1 / (2 * np.sin(x))), np.pi / 6, 5 * np.pi / 6)
def a3w():
    # M.re: asin of an argument rounded past 1 at a node next to the end is
    # complex (Im ~ 1e-140, a rounding artefact); the area is its real part.
    return M.re(M.quad(lambda x: M.pi - 2 * M.asin(1 / (2 * M.sin(x))), [M.pi / 6, 5 * M.pi / 6]))
# 4
def a4m():
    G = lambda y: -np.log(0.5 + y * y / 8) - y * y
    s = lambda y: np.sqrt(-np.log(0.5 + y * y / 8))
    yA = nt(G, -0.75); yB = nt(G, 0.75)
    return (q(lambda y: 2 * np.sqrt(y + s(y)), yA, yB, [0.0])
            + q(lambda y: 2 * (np.sqrt(y + s(y)) - np.sqrt(y - s(y))), yB, 2))
def a4w():
    g = lambda y: M.mpf(1) / 2 + y ** 2 / 8
    s = lambda y: M.sqrt(-M.log(g(y)))
    yA = M.findroot(lambda y: -M.log(g(y)) - y ** 2, mf("-0.75")); yB = M.findroot(lambda y: -M.log(g(y)) - y ** 2, mf("0.75"))
    return (M.quad(lambda y: 2 * M.sqrt(y + s(y)), [yA, 0, yB])
            + M.quad(lambda y: 2 * (M.sqrt(y + s(y)) - M.sqrt(y - s(y))), [yB, 2]))
# 5
def a5m():
    p = lambda x: x ** 3 - 3 * x + 1
    r1 = nt(p, -1.9); r2 = nt(p, 0.35)
    return q(lambda x: 2 * np.sqrt(p(x)), r1, r2)
def a5w():
    p = lambda x: x ** 3 - 3 * x + 1
    r1 = M.findroot(p, mf("-1.9")); r2 = M.findroot(p, mf("0.35"))
    return M.quad(lambda x: 2 * M.sqrt(p(x)), [r1, r2])
# 6
def a6m():
    x0 = nt(lambda x: np.cos(x) - x * x, 0.8)
    return 4 * q(lambda x: (np.cos(x) - x * x) ** 0.25, 0, x0)
def a6w():
    x0 = M.findroot(lambda x: M.cos(x) - x ** 2, mf("0.8"))
    return 4 * M.quad(lambda x: M.root(M.cos(x) - x ** 2, 4), [0, x0])
# 7
def a7m():
    h = lambda u: 1 + np.sin(np.sqrt(2) * u) - u * u
    u1 = nt(h, -0.55); u2 = nt(h, 1.4)
    return q(lambda u: 2 * np.sqrt(h(u)), u1, u2)
def a7w():
    h = lambda u: 1 + M.sin(M.sqrt(2) * u) - u ** 2
    u1 = M.findroot(h, mf("-0.55")); u2 = M.findroot(h, mf("1.4"))
    return M.quad(lambda u: 2 * M.sqrt(h(u)), [u1, u2])
# 8
def a8m():
    g = lambda y: y * y / 16 + np.exp(np.cos(y) - 1)
    w = lambda y: np.cbrt(np.sqrt(-np.log(g(y))) - 2 * y) - np.cbrt(-np.sqrt(-np.log(g(y))) - 2 * y)
    m = nt(lambda y: g(y) - 1, -3.7); p = nt(lambda y: g(y) - 1, 3.7)
    return q(w, m, p, [0.0])
def a8w():
    g = lambda y: y ** 2 / 16 + M.exp(M.cos(y) - 1)
    w = lambda y: rc(M.sqrt(-M.log(g(y))) - 2 * y) - rc(-M.sqrt(-M.log(g(y))) - 2 * y)
    m = M.findroot(lambda y: g(y) - 1, mf("-3.7")); p = M.findroot(lambda y: g(y) - 1, mf("3.7"))
    return M.quad(w, [m, 0, p])
# 9
def a9m():
    c = lambda x: np.sqrt(2.25 - x * x); a = lambda x: np.arccos(1 - np.cos(x))
    x1 = nt(lambda x: c(x) - a(x), 0.75); x2 = nt(lambda x: c(x) - a(x), 1.3)
    return 4 * (q(c, 0, x1) + q(a, x1, x2) + q(c, x2, 1.5))
def a9w():
    c = lambda x: M.sqrt(M.mpf(9) / 4 - x ** 2); a = lambda x: M.acos(1 - M.cos(x))
    x1 = M.findroot(lambda x: c(x) - a(x), mf("0.75")); x2 = M.findroot(lambda x: c(x) - a(x), mf("1.3"))
    # M.re: findroot's secant steps leave x2 with an Im ~ 1e-100 artefact
    return M.re(4 * (M.quad(c, [0, x1]) + M.quad(a, [x1, x2]) + M.quad(c, [x2, M.mpf(3) / 2])))
# 10
S10 = ["-7.94", "-7.04", "-3.04", "3.25", "5.85"]
def a10m():
    g = lambda y: y * y / 100 + np.exp(np.sin(y))
    w = lambda y: np.cbrt(-y + np.sqrt(-np.log(g(y)))) - np.cbrt(-y - np.sqrt(-np.log(g(y))))
    y1, y2, y3, y4, y5 = [nt(lambda y: g(y) - 1, float(s)) for s in S10]
    yc = nt(lambda y: -np.log(g(y)) - y * y, -0.86)
    return q(w, y1, y2) + q(w, y3, 0, [yc]) + q(w, y4, y5)
def a10w():
    g = lambda y: y ** 2 / 100 + M.exp(M.sin(y))
    w = lambda y: rc(-y + M.sqrt(-M.log(g(y)))) - rc(-y - M.sqrt(-M.log(g(y))))
    y1, y2, y3, y4, y5 = [M.findroot(lambda y: g(y) - 1, mf(s)) for s in S10]
    yc = M.findroot(lambda y: -M.log(g(y)) - y ** 2, mf("-0.86"))
    return M.quad(w, [y1, y2]) + M.quad(w, [y3, yc, 0]) + M.quad(w, [y4, y5])

for k, (fm, fw) in enumerate([(a1m, a1w), (a2m, a2w), (a3m, a3w), (a4m, a4w), (a5m, a5w),
                              (a6m, a6w), (a7m, a7w), (a8m, a8w), (a9m, a9w), (a10m, a10w)], 1):
    vrun("bd7.%d machine" % k, fm, 10, 3, cap=60)
    vrun("bd7.%d wp50" % k, fw, 48, cap=60, dps=60)
