"""Experiment 226 -- SIAM 100-digit holdout B for Problem 1 (Python column).

The integrals of oscillatory_tails.m as written.  Machine: scipy.integrate.quad
(limit=1000); where the integrand is (non-oscillatory) x (sin|cos)(x) on
[a, oo) -- bc1.5, bc1.7, bc1.10 -- the idiomatic scipy spelling is quad's
Fourier weight (weight='sin'/'cos', wvar=1: QUADPACK QAWF).  wp50 (dps 60):
mpmath.quad on [0, 1] for the x -> 0 problems; mpmath.quadosc on [a, oo)
(omega=1 for a period-2pi oscillation, else the documented zeros= form)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath as mp
from scipy.integrate import quad
from scipy.special import j1

require(["scipy.integrate:quad", "mpmath:quad", "mpmath:quadosc"])
inf = float("inf")


def pair(k, mach, wp):
    vrun("bc1.%d machine" % k, mach, 10, 3, cap=60)
    vrun("bc1.%d wp50" % k, wp, 48, cap=120, dps=60)


pair(1, lambda: quad(lambda x: math.sin(1 / x), 0, 1, limit=1000)[0],
     lambda: mp.quad(lambda x: mp.sin(1 / x), [0, 1]))
pair(2, lambda: quad(lambda x: math.sin(1 / x) / x, 0, 1, limit=1000)[0],
     lambda: mp.quad(lambda x: mp.sin(1 / x) / x, [0, 1]))
pair(3, lambda: quad(lambda x: math.cos(math.log(x) / x), 0, 1, limit=1000)[0],
     lambda: mp.quad(lambda x: mp.cos(mp.log(x) / x), [0, 1]))
pair(4, lambda: quad(lambda x: math.sin(x * x), 0, inf, limit=1000)[0],
     lambda: mp.quadosc(lambda x: mp.sin(x**2), [0, mp.inf], zeros=lambda n: mp.sqrt(n * mp.pi)))
pair(5, lambda: quad(lambda x: 1 / math.sqrt(x), 0, inf, weight="cos", wvar=1)[0],
     lambda: mp.quadosc(lambda x: mp.cos(x) / mp.sqrt(x), [0, mp.inf], omega=1))
pair(6, lambda: quad(lambda x: math.sin(1 / x) / math.sqrt(x), 0, 1, limit=1000)[0],
     lambda: mp.quad(lambda x: mp.sin(1 / x) / mp.sqrt(x), [0, 1]))
pair(7, lambda: quad(lambda x: 1 / (x + 1), 0, inf, weight="sin", wvar=1)[0],
     lambda: mp.quadosc(lambda x: mp.sin(x) / (x + 1), [0, mp.inf], omega=1))
pair(8, lambda: quad(lambda x: j1(x) / math.sqrt(x), 0, inf, limit=1000)[0],
     lambda: mp.quadosc(lambda x: mp.besselj(1, x) / mp.sqrt(x), [0, mp.inf], omega=1))
pair(9, lambda: quad(lambda x: math.cos(x**3), 0, inf, limit=1000)[0],
     lambda: mp.quadosc(lambda x: mp.cos(x**3), [0, mp.inf],
                        zeros=lambda n: mp.cbrt((n - mp.mpf(1) / 2) * mp.pi)))
pair(10, lambda: quad(lambda x: 1 / math.log(x), 2, inf, weight="sin", wvar=1)[0],
     lambda: mp.quadosc(lambda x: mp.sin(x) / mp.log(x), [2, mp.inf], omega=1))
