"""Experiment 239 -- SIAM 100-Digit Challenge, holdout B for Appendix D
problem 4 (Python column).  Ten roots of special functions / series; see
series_roots.m.  bd4.1 machine is scipy (brentq on scipy.special.expi(log x),
li(x) = Ei(log x)); every other row is mpmath.findroot (secant, same start as
the .m file) on the same formulation -- the series rows via mpmath.nsum, mu(k)
from sympy.mobius, the built-ins mpmath.riemannr / siegelz / primezeta -- at
mp.dps 15 for the machine rows and 60 for wp50.  The answer for bd4.4/5/9 is
t = log x."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import mpmath
from mpmath import mpf, inf
from scipy.optimize import brentq
from scipy.special import expi
from sympy import mobius
import math

require(["scipy.optimize:brentq", "scipy.special:expi", "mpmath:findroot", "mpmath:nsum",
         "mpmath:riemannr", "mpmath:siegelz", "mpmath:primezeta", "sympy:mobius"])

fr = mpmath.findroot
ns = mpmath.nsum
cases = {
    1: lambda: fr(mpmath.li, mpf(3) / 2),
    2: lambda: fr(lambda s: ns(lambda n: n ** -s, [1, inf]) - 2, mpf(17) / 10),
    3: lambda: fr(lambda s: ns(lambda n: n ** -s / mpmath.log(n), [2, inf]) - 1, mpf(8) / 5),
    4: lambda: fr(lambda t: mpmath.riemannr(mpmath.exp(t)) - mpf(10) ** -3, -28),
    5: lambda: fr(lambda t: ns(lambda k: mobius(int(k)) / k * mpmath.li(mpmath.exp(t / k)),
                               [1, inf]) - mpf(1) / 10, -5),
    6: lambda: fr(lambda x: ns(lambda k: mobius(int(k)) * x ** k, [1, inf]), mpf(3) / 5),
    7: lambda: fr(mpmath.siegelz, 14),
    8: lambda: fr(lambda s: mpmath.primezeta(s) - 1, mpf(7) / 5),
    9: lambda: fr(lambda t: mpmath.riemannr(mpmath.exp(t)) - mpf(10) ** -9, -2800),
    10: lambda: fr(lambda s: ns(lambda n: -mpmath.log(n) * n ** -s, [1, inf]) + 1, mpf(19) / 10),
}
cap = {2: 120, 3: 120, 5: 120, 6: 120, 9: 120, 10: 120}
vrun("bd4.1 machine", lambda: brentq(lambda x: expi(math.log(x)), 1.2, 2.0, xtol=1e-15, rtol=8.9e-16),
     10, cap=60)
vrun("bd4.1 wp50", cases[1], 48, cap=60, dps=60)
for p in range(2, 11):
    vrun("bd4.%d machine" % p, lambda: float(cases[p]()), 10, cap=60, dps=15)
    vrun("bd4.%d wp50" % p, cases[p], 48, cap=cap.get(p, 60), dps=60)
