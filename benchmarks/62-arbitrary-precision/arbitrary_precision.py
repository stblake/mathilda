#!/usr/bin/env python3
"""Experiment 62 -- Arbitrary-precision numerics (mpmath column).

Same kernels as ``arbitrary_precision.m``, same order and labels.

Baseline is mpmath -- the standard Python arbitrary-precision library (scipy is
double-only, so it cannot appear here).  mpmath is a weaker baseline than scipy
in one sense: pure Python unless gmpy2 is installed, in which case its core is
GMP/MPFR -- the same library Mathilda's precision path uses.  The .py docstring
of a case says which.

Single-run (bench_once): these computations are expensive on their own.  Checks
use a deep 6-digit window win(v, s) = digits s-5..s of v, matching the .m side.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import mpmath as mp

from harness import bench_once, check, require

require(["mpmath:pi", "mpmath:nsum", "mpmath:quad", "mpmath:findroot"])


def win(v, start):
    return int(mp.floor(v * mp.mpf(10) ** start)) % 10 ** 6


def n_pi():
    mp.mp.dps = 10020
    return +mp.pi

bench_once("N[Pi, 10000]", n_pi)
check("N[Pi, 10000]", win(n_pi(), 5000))


def n_gamma():
    mp.mp.dps = 1020
    return mp.gamma(mp.mpf(1) / 3)

bench_once("N[Gamma[1/3], 1000]", n_gamma)
check("N[Gamma[1/3], 1000]", win(n_gamma(), 500))


def n_log2():
    mp.mp.dps = 5020
    return mp.log(2)

bench_once("N[Log[2], 5000]", n_log2)
check("N[Log[2], 5000]", win(n_log2(), 2000))


def s_zeta2():
    mp.mp.dps = 120
    return mp.nsum(lambda k: 1 / k ** 2, [1, mp.inf])

bench_once("NSum 1/k^2 WorkingPrecision 100", s_zeta2)
check("NSum 1/k^2 WorkingPrecision 100", win(s_zeta2(), 60))


def i_gauss():
    mp.mp.dps = 220
    return mp.quad(lambda x: mp.e ** (-x * x), [0, 1])

bench_once("NIntegrate Exp[-x^2] WorkingPrecision 200", i_gauss)
check("NIntegrate Exp[-x^2] WorkingPrecision 200", win(i_gauss(), 50))


def r_dottie():
    mp.mp.dps = 220
    return mp.findroot(lambda x: mp.cos(x) - x, 1)

bench_once("FindRoot Cos[x]=x WorkingPrecision 200", r_dottie)
check("FindRoot Cos[x]=x WorkingPrecision 200", win(r_dottie(), 50))
