"""Experiment 241 -- holdout B for SIAM Appendix D Problem 6 (Python column).

The ten integrals of oscillatory_towers.m as written: scipy.integrate.quad at
machine precision (tight tolerances, limit 200), mpmath.quad (tanh-sinh, no
splitting) at 60 digits.  See oscillatory_towers.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath
from scipy.integrate import quad

require(["scipy.integrate:quad", "mpmath:quad"])

M = mpmath
cases = [
  (lambda x: math.exp(math.sin(math.exp(x))), lambda x: M.exp(M.sin(M.exp(x))), 0, 1),
  (lambda x: math.cos(x + math.exp(math.exp(x + 0.5))), lambda x: M.cos(x + M.exp(M.exp(x + M.mpf(1)/2))), -1, 1),
  (lambda x: math.exp(math.cos(x**3 + math.exp(x))), lambda x: M.exp(M.cos(x**3 + M.exp(x))), 0, 2),
  (lambda x: math.log(2 + math.sin(math.exp(math.exp(x)))), lambda x: M.log(2 + M.sin(M.exp(M.exp(x)))), -1, 1),
  (lambda x: 1 / (1.25 + math.cos(math.exp(4 * x))), lambda x: 1 / (M.mpf(5)/4 + M.cos(M.exp(4 * x))), 0, 1),
  (lambda x: math.sin(math.exp(math.exp(2 * x))), lambda x: M.sin(M.exp(M.exp(2 * x))), 0, 1),
  (lambda x: math.cos(50 * math.sin(x) - 3 * x), lambda x: M.cos(50 * M.sin(x) - 3 * x), 0, "pi"),
  (lambda x: math.exp(x) * math.cos(math.sinh(math.exp(x + 0.5))), lambda x: M.exp(x) * M.cos(M.sinh(M.exp(x + M.mpf(1)/2))), -1, 1),
  (lambda x: math.cos(math.exp(math.exp(math.exp(x)))), lambda x: M.cos(M.exp(M.exp(M.exp(x)))), -1, "1/2"),
  (lambda x: math.exp(x + math.sin(math.exp(math.exp(x + 0.5)))), lambda x: M.exp(x + M.sin(M.exp(M.exp(x + M.mpf(1)/2)))), -1, 1),
]
def endpoint_f(b):
    return {"pi": math.pi, "1/2": 0.5}.get(b, b)
def endpoint_m(b):
    if b == "pi":
        return +M.pi
    if b == "1/2":
        return M.mpf(1) / 2
    return M.mpf(b)

for k, (fm, fw, a, b) in enumerate(cases, 1):
    vrun("bd6.%d machine" % k,
         lambda: quad(fm, a, endpoint_f(b), epsabs=1e-14, epsrel=1e-13, limit=200)[0], 10, 3, cap=60)
    vrun("bd6.%d wp50" % k,
         lambda: M.quad(fw, [M.mpf(a), endpoint_m(b)]), 48, cap=60, dps=60)
