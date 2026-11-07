"""Experiment 209 -- SIAM 100-digit holdout A, family d4 (Python column).

scipy has neither li-roots of this kind nor a series summer, so every row is
mpmath.findroot on the same formulation as number_theory_roots_family.m:
built-ins mpmath.li / riemannr / siegelz where the .m file uses LogIntegral /
RiemannR / RiemannSiegelZ, and mpmath.nsum where it uses NSum (mu(k) from
sympy.mobius, Lambda(n) from mpmath.mangoldt).  dps 15 for the machine rows,
dps 60 for wp50.  Per-case cap 60 s (120 s for hd4.10 wp50)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import mpmath
from mpmath import mpf
from sympy import mobius

require(["mpmath:findroot", "mpmath:nsum", "mpmath:li", "mpmath:riemannr", "mpmath:siegelz",
         "mpmath:mangoldt", "sympy:mobius"])
fr, nsum, inf = mpmath.findroot, mpmath.nsum, mpmath.inf


def gram(x):
    L = mpmath.log(x)
    return 1 + nsum(lambda k: L ** k / (k * mpmath.factorial(k) * mpmath.zeta(k + 1)), [1, inf])


def moeb(x):
    return nsum(lambda k: mobius(int(k)) / k * mpmath.li(x ** (1 / k)), [1, inf])


def imu(s):
    return nsum(lambda k: mobius(int(k)) / k ** s, [1, inf])


def lam(s):
    return nsum(lambda k: mpmath.mangoldt(int(k)) / k ** s, [1, inf])


CASES = {
    1: lambda: fr(mpmath.li, 2),
    2: lambda: fr(lambda x: mpmath.riemannr(x) - 100, 500),
    3: lambda: fr(lambda x: gram(x) - 2, 3),
    4: lambda: fr(lambda x: moeb(x) - 1000, 8000),
    5: lambda: fr(lambda x: gram(x) - mpf(1) / 2, mpf(1) / 4),
    6: lambda: fr(lambda s: imu(s) - mpf(1) / 2, 2),
    7: lambda: fr(lambda s: lam(s) - 1, 2),
    8: lambda: fr(mpmath.siegelz, 14),
    9: lambda: fr(lambda x: gram(x) - mpf(10) ** 10, 2e11),
    10: lambda: mpmath.exp(fr(lambda t: mpmath.riemannr(mpmath.exp(t)), -35300)),
}
for k in range(1, 11):
    vrun("hd4.%d machine" % k, CASES[k], 10, cap=60, dps=15)
    vrun("hd4.%d wp50" % k, CASES[k], 48, cap=120 if k == 10 else 60, dps=60)
