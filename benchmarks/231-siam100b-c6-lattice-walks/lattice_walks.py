"""Experiment 231 -- SIAM 100-digit HOLDOUT B for Problem 6 (Python column).

The same ten statements as lattice_walks.m: each lattice Green's function as
the folded Fourier integral over [0,pi]^2 (or [0,pi]^3) and, where asked, a
root solve on it.  Machine rows: scipy.integrate.dblquad / tplquad (epsabs 0,
epsrel 1e-12) inside scipy.optimize.brentq on the same bracket the .m half
uses.  wp50 rows: mpmath.quad (2-D / 3-D tanh-sinh) inside mpmath.findroot
(secant from the bracket ends), dps 60.  See lattice_walks.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import warnings
warnings.simplefilter("ignore")      # dblquad's IntegrationWarning
import mpmath
from scipy.integrate import dblquad, tplquad
from scipy.optimize import brentq

require(["scipy.integrate:dblquad", "scipy.integrate:tplquad", "scipy.optimize:brentq",
         "mpmath:quad", "mpmath:findroot"])

PI = math.pi
MC, WC = 90, 75          # per-case caps (s): machine, wp50


def re2(c, s):
    return (1 - c) / ((1 - c) ** 2 + s * s)


# integrands f(k1, k2[, k3], e, M) with M = math or mpmath, and their normalisation
def f1(k1, k2, e, M): return re2((M.cos(k1) + M.cos(k2)) / 2, 2 * e * M.sin(k1))
def f2(k1, k2, e, M):
    c = (M.cos(k1) + M.cos(k2)) / 2
    return re2(c, 2 * e * (M.sin(k1) + M.sin(k2))) + re2(c, 2 * e * (M.sin(k1) - M.sin(k2)))
def f3(k1, k2, e, M): return re2((2 * M.cos(k1) + M.cos(k2)) / 3, 2 * e * M.sin(k1))
def f4(k1, k2, e, M):
    return (re2((M.cos(k1) + M.cos(k2) + M.cos(k1 + k2)) / 3, 2 * e * M.sin(k1))
            + re2((M.cos(k1) + M.cos(k2) + M.cos(k1 - k2)) / 3, 2 * e * M.sin(k1)))
def f6(k1, k2, e, M): return re2((M.cos(k1) + M.cos(k2) + 2 * M.cos(k1) * M.cos(k2)) / 4, 2 * e * M.sin(k1))
def h9(k1, k2, e, M):
    c = (M.cos(k1) + M.cos(k2)) / 2; s = 2 * e * M.sin(k1)
    return (M.cos(k1) * (1 - c) + s * M.sin(k1)) / ((1 - c) ** 2 + s * s)
def f10(k1, k2, e, M):
    c = M.cos(k1) * M.cos(k2)
    return re2(c, 2 * e * M.sin(k1 + k2)) + re2(c, 2 * e * M.sin(k1 - k2))
def f5(k1, k2, k3, e, M): return 1 / (1 - (M.cos(k1) + M.cos(k2) + M.cos(k3)) / 3)
def f8(k1, k2, k3, e, M): return re2((M.cos(k1) + M.cos(k2) + M.cos(k3)) / 3, 2 * e * M.sin(k1))


def Gd(f, e, norm):        # machine 2-D
    return dblquad(lambda k2, k1: f(k1, k2, e, math), 0, PI, 0, PI, epsabs=0, epsrel=1e-12)[0] / norm


def Gt(f, e):              # machine 3-D
    return tplquad(lambda k3, k2, k1: f(k1, k2, k3, e, math), 0, PI, 0, PI, 0, PI,
                   epsabs=0, epsrel=1e-12)[0] / PI ** 3


def Gm(f, e, norm):        # mpmath 2-D
    p = mpmath.pi
    return mpmath.quad(lambda k1, k2: f(k1, k2, e, mpmath), [0, p], [0, p]) / norm(p)


def Gm3(f, e):
    p = mpmath.pi
    return mpmath.quad(lambda k1, k2, k3: f(k1, k2, k3, e, mpmath), [0, p], [0, p], [0, p]) / p ** 3


one = lambda p: p ** 2
two = lambda p: 2 * p ** 2
mpf = mpmath.mpf


def mroot(f, norm, target, lo, hi):
    return lambda: brentq(lambda e: Gd(f, e, norm(PI)) - target, lo, hi, xtol=1e-16, rtol=1e-15)


def wroot(f, norm, target, lo, hi):
    return lambda: mpmath.findroot(lambda e: Gm(f, e, norm) - target, (mpf(lo), mpf(hi)))


vrun("bc6.1 machine", mroot(f1, one, 2.5, 0.02, 0.04), 10, cap=MC)
vrun("bc6.1 wp50", wroot(f1, one, mpf(5) / 2, "0.02", "0.04"), 48, cap=WC, dps=60)
vrun("bc6.2 machine", mroot(f2, two, 2, 0.03, 0.06), 10, cap=MC)
vrun("bc6.2 wp50", wroot(f2, two, 2, "0.03", "0.06"), 48, cap=WC, dps=60)
vrun("bc6.3 machine", mroot(f3, one, 2, 0.06, 0.1), 10, cap=MC)
vrun("bc6.3 wp50", wroot(f3, one, 2, "0.06", "0.1"), 48, cap=WC, dps=60)
vrun("bc6.4 machine", mroot(f4, two, 2, 0.03, 0.06), 10, cap=MC)
vrun("bc6.4 wp50", wroot(f4, two, 2, "0.03", "0.06"), 48, cap=WC, dps=60)
vrun("bc6.5 machine", lambda: 1 - 1 / Gt(f5, 0), 10, cap=MC)
vrun("bc6.5 wp50", lambda: 1 - 1 / Gm3(f5, 0), 48, cap=WC, dps=60)
vrun("bc6.6 machine", mroot(f6, one, 2, 0.02, 0.05), 10, cap=MC)
vrun("bc6.6 wp50", wroot(f6, one, 2, "0.02", "0.05"), 48, cap=WC, dps=60)
vrun("bc6.7 machine", lambda: 1 - 1 / Gd(f1, 0.1, PI ** 2), 10, 3, cap=MC)
vrun("bc6.7 wp50", lambda: 1 - 1 / Gm(f1, mpf(1) / 10, one), 48, cap=WC, dps=60)
vrun("bc6.8 machine", lambda: brentq(lambda e: Gt(f8, e) - 4 / 3, 0.05, 0.1, xtol=1e-16, rtol=1e-15), 10, cap=MC)
vrun("bc6.8 wp50", lambda: mpmath.findroot(lambda e: Gm3(f8, e) - mpf(4) / 3, (mpf("0.05"), mpf("0.1"))),
     48, cap=WC, dps=60)
vrun("bc6.9 machine", lambda: Gd(h9, 0.05, PI ** 2) / Gd(f1, 0.05, PI ** 2), 10, 3, cap=MC)
vrun("bc6.9 wp50", lambda: Gm(h9, mpf(1) / 20, one) / Gm(f1, mpf(1) / 20, one), 48, cap=WC, dps=60)
vrun("bc6.10 machine", mroot(f10, two, 5 / 3, 0.08, 0.13), 10, cap=MC)
vrun("bc6.10 wp50", wroot(f10, two, mpf(5) / 3, "0.08", "0.13"), 48, cap=WC, dps=60)
