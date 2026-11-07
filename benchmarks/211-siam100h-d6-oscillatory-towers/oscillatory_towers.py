"""Experiment 211 -- SIAM 100-digit HOLDOUT A, family D6 (Python column).

The ten integrals of oscillatory_towers.m, as written: scipy.integrate.quad
at machine precision (epsabs 1e-14, epsrel 1e-13, limit 1000), mpmath.quad
(tanh-sinh, no splitting) at 60 digits."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath as M
from scipy.integrate import quad

require(["scipy.integrate:quad", "mpmath:quad"])

m, e, pi = math, math.exp, math.pi
F = {  # k: (machine integrand, mpmath integrand, a, b)
  1: (lambda x: m.cos(e(x*x)),                 lambda x: M.cos(M.exp(x*x)), 0, 2),
  2: (lambda x: e(x + m.cos(e(e(x)))),         lambda x: M.exp(x + M.cos(M.exp(M.exp(x)))), -1, 1),
  3: (lambda x: m.sin(e(e(2*x))),              lambda x: M.sin(M.exp(M.exp(2*x))), 0, 1),
  4: (lambda x: m.sin(x**3),                   lambda x: M.sin(x**3), 0, 3),
  5: (lambda x: e(m.sin(e(x))),                lambda x: M.exp(M.sin(M.exp(x))), 0, 'pi'),
  6: (lambda x: m.sin(e(4*x))/(2 + x),         lambda x: M.sin(M.exp(4*x))/(2 + x), -1, 1),
  7: (lambda x: x*e(x*x)*m.sin(e(x*x)),        lambda x: x*M.exp(x*x)*M.sin(M.exp(x*x)), 0, 2),
  8: (lambda x: e(x + m.sin(e(e(x + 0.5)))),   lambda x: M.exp(x + M.sin(M.exp(M.exp(x + M.mpf(1)/2)))), -1, 1),
  9: (lambda x: m.cos(x*x + x),                lambda x: M.cos(x*x + x), 0, 8),
  10: (lambda x: m.cos(1/(x + 0.01)),          lambda x: M.cos(1/(x + M.mpf(1)/100)), 0, 1),
}
for k in range(1, 11):
    f, g, a, b = F[k]
    bm = pi if b == 'pi' else b
    vrun("hd6.%d machine" % k,
         lambda: quad(f, a, bm, epsabs=1e-14, epsrel=1e-13, limit=1000)[0], 10, 3, cap=30)
    vrun("hd6.%d wp50" % k,
         lambda: M.quad(g, [a, M.pi if b == 'pi' else b]), 48, cap=60, dps=60)
