"""Experiment 201 -- SIAM 100-Digit Challenge, holdout A for Problem 6 (Python column).

The same Brillouin-zone integrals as lattice_return.m: scipy.integrate
dblquad / tplquad (epsabs 0, epsrel 1e-12) inside
scipy.optimize.brentq at machine precision; mpmath.quad (tanh-sinh, 1-, 2- or
3-D) inside mpmath.findroot (secant from the .m file's two starts) at dps 60.
Every case is capped at 60 s.  See lattice_return.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import warnings
warnings.simplefilter("ignore")
import mpmath
from scipy.integrate import quad, dblquad, tplquad
from scipy.optimize import brentq

require(["scipy.integrate:dblquad", "scipy.integrate:tplquad", "scipy.optimize:brentq",
         "mpmath:quad", "mpmath:findroot"])

PI = math.pi
R = lambda c, s: (1 - c) / ((1 - c) ** 2 + s * s)


# ---- machine: integrands of (k1, k2) -> value; dblquad wants f(y, x)
def dq(f, sym):
    lo2 = 0 if sym else -PI
    v = dblquad(lambda k2, k1: f(k1, k2), 0, PI, lo2, PI, epsabs=0, epsrel=1e-12)[0]
    return v / PI ** 2 if sym else v / (2 * PI ** 2)


cS = lambda k1, k2: (math.cos(k1) + math.cos(k2)) / 2
gEW = lambda e: dq(lambda k1, k2: R(cS(k1, k2), 2 * e * math.sin(k1)), True)
gDiag = lambda e: dq(lambda k1, k2: R(cS(k1, k2), 2 * e * (math.sin(k1) + math.sin(k2))), False)
gAni = lambda e: dq(lambda k1, k2: R(2 * math.cos(k1) / 3 + math.cos(k2) / 3, 2 * e * math.sin(k1)), True)
gEW3 = lambda e: tplquad(lambda k3, k2, k1: R((math.cos(k1) + math.cos(k2) + math.cos(k3)) / 3, 2 * e * math.sin(k1)),
                         0, PI, 0, PI, 0, PI, epsabs=0, epsrel=1e-12)[0] / PI ** 3
gTri = lambda e: dq(lambda k1, k2: R((math.cos(k1) + math.cos(k2) + math.cos(k1 + k2)) / 3, 2 * e * math.sin(k1)), False)
g9 = lambda: dq(lambda k1, k2: R(cS(k1, k2), math.sin(k1) / 10 + math.sin(k2) / 5), False)
u3 = lambda: tplquad(lambda k3, k2, k1: 1 / (1 - (math.cos(k1) + math.cos(k2) + math.cos(k3)) / 3),
                     0, PI, 0, PI, 0, PI, epsabs=0, epsrel=1e-12)[0] / PI ** 3

# ---- mpmath
mc, ms, mpi = mpmath.cos, mpmath.sin, lambda: mpmath.pi


def mq(f, sym):
    lo2 = 0 if sym else -mpi()
    v = mpmath.quad(f, [0, mpi()], [lo2, mpi()])
    return v / mpi() ** 2 if sym else v / (2 * mpi() ** 2)


mcS = lambda k1, k2: (mc(k1) + mc(k2)) / 2
GEW = lambda e: mq(lambda k1, k2: R(mcS(k1, k2), 2 * e * ms(k1)), True)
GDiag = lambda e: mq(lambda k1, k2: R(mcS(k1, k2), 2 * e * (ms(k1) + ms(k2))), False)
GAni = lambda e: mq(lambda k1, k2: R(2 * mc(k1) / 3 + mc(k2) / 3, 2 * e * ms(k1)), True)
GEW3 = lambda e: mpmath.quad(lambda k1, k2, k3: R((mc(k1) + mc(k2) + mc(k3)) / 3, 2 * e * ms(k1)),
                             [0, mpi()], [0, mpi()], [0, mpi()]) / mpi() ** 3
GTri = lambda e: mq(lambda k1, k2: R((mc(k1) + mc(k2) + mc(k1 + k2)) / 3, 2 * e * ms(k1)), False)
G9 = lambda: mq(lambda k1, k2: R(mcS(k1, k2), ms(k1) / 10 + ms(k2) / 5), False)
U3 = lambda: mpmath.quad(lambda k1, k2, k3: 1 / (1 - (mc(k1) + mc(k2) + mc(k3)) / 3),
                         [0, mpi()], [0, mpi()], [0, mpi()]) / mpi() ** 3

mf = mpmath.mpf
# (machine G, mpmath G, target G, bracket)
roots = {2: (gEW, GEW, 5 / 3, mf(5) / 3, (0.10, 0.11)),
         3: (gEW, GEW, 5 / 2, mf(5) / 2, (0.027, 0.029)),
         4: (gDiag, GDiag, 2, 2, (0.04, 0.05)),
         5: (gAni, GAni, 2, 2, (0.075, 0.085)),
         6: (gEW, GEW, 5 / 4, mf(5) / 4, (0.2, 0.22)),
         8: (gEW3, GEW3, 10 / 7, mf(10) / 7, (0.03, 0.035)),
         10: (gTri, GTri, 2, 2, (0.045, 0.05))}


def brk(f, T, a, b):
    """brentq on the bracket (the .m file's two FindRoot starts)."""
    return brentq(lambda e: f(e) - T, a, b, xtol=1e-16, rtol=1e-15)


for k in range(1, 11):
    lm, lw = "hc6.%d machine" % k, "hc6.%d wp50" % k
    if k == 1:
        vrun(lm, lambda: 1 - 1 / gEW(0.1), 10, cap=60)
        vrun(lw, lambda: 1 - 1 / GEW(mf(1) / 10), 48, cap=60, dps=60)
    elif k == 7:
        vrun(lm, lambda: 1 - 1 / u3(), 10, cap=60)
        vrun(lw, lambda: 1 - 1 / U3(), 48, cap=60, dps=60)
    elif k == 9:
        vrun(lm, lambda: 1 - 1 / g9(), 10, cap=60)
        vrun(lw, lambda: 1 - 1 / G9(), 48, cap=60, dps=60)
    else:
        f, F, T, TW, (a, b) = roots[k]
        vrun(lm, lambda f=f, T=T, a=a, b=b: brk(f, T, a, b), 10, cap=60)
        vrun(lw, lambda F=F, TW=TW, a=a, b=b: mpmath.findroot(
            lambda e: F(e) - TW,
            (mf(a), mf(b))), 48, cap=60, dps=60)
