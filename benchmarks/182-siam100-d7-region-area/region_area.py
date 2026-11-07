"""Experiment 182 -- SIAM 100-Digit Challenge, Appendix D Problem 7 (Python column).

The 1-D statement of region_area.m: island ends and the interior kink by
root-finding from the same starting values (scipy.optimize.newton, secant
mode; mpmath.findroot), then the x-interval length w(y) integrated over
[y1, yc, 0] and [y2, y3] (scipy.integrate.quad with the kink as a
breakpoint; mpmath.quad on the same subintervals at 60 digits).
numpy.cbrt is the real cube root; mpmath.cbrt is the principal complex one,
so the wp50 row takes the real cube root as sign(t) * cbrt(|t|)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import numpy as np
import mpmath
from scipy.optimize import newton
from scipy.integrate import quad

require(["scipy.optimize:newton", "scipy.integrate:quad", "mpmath:findroot", "mpmath:quad"])

g = lambda y: y * y / 32 + np.exp(np.sin(y))
def w(y):
    s = np.sqrt(-np.log(g(y)))
    return np.cbrt(-y + s) - np.cbrt(-y - s)
def area_machine():
    y1 = newton(lambda y: g(y) - 1, -2.8, tol=1e-15)
    yc = newton(lambda y: -np.log(g(y)) - y * y, -0.8, tol=1e-15)
    y2 = newton(lambda y: g(y) - 1, 3.8, tol=1e-15)
    y3 = newton(lambda y: g(y) - 1, 4.45, tol=1e-15)
    return (quad(w, y1, 0, points=[yc], epsabs=1e-14, epsrel=1e-13, limit=200)[0]
            + quad(w, y2, y3, epsabs=1e-14, epsrel=1e-13, limit=200)[0])
vrun("d7.machine", area_machine, 10, 3)

G = lambda y: y ** 2 / 32 + mpmath.exp(mpmath.sin(y))
rc = lambda t: mpmath.sign(t) * mpmath.cbrt(abs(t))
def W(y):
    s = mpmath.sqrt(-mpmath.log(G(y)))
    return rc(-y + s) - rc(-y - s)
def area_wp50():
    y1 = mpmath.findroot(lambda y: G(y) - 1, mpmath.mpf("-2.8"))
    yc = mpmath.findroot(lambda y: -mpmath.log(G(y)) - y ** 2, mpmath.mpf("-0.8"))
    y2 = mpmath.findroot(lambda y: G(y) - 1, mpmath.mpf("3.8"))
    y3 = mpmath.findroot(lambda y: G(y) - 1, mpmath.mpf("4.45"))
    return mpmath.quad(W, [y1, yc, 0]) + mpmath.quad(W, [y2, y3])
vrun("d7.wp50", area_wp50, 48, dps=60)
