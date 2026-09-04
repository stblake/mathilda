#!/usr/bin/env python3
"""Experiment 55 -- Special functions over machine arrays (scipy.special column).

Same kernels as ``vectorized_special_functions.m``, same order and labels.

EXECUTION comparison: scipy.special is compiled C (SIMD ufuncs).  For the
vectorized Mathilda heads the gap is overhead; for BesselI/BesselK -- which have
no vector kernel in Mathilda and thread scalar-by-scalar -- the gap is a missing
kernel, and the ratio is large.  All ten are one ufunc call here.

Normalization: Mathematica's FresnelC/FresnelS and scipy.special.fresnel both
use the pi*t^2/2 argument; sici returns (Si, Ci); fresnel returns (S, C).
Checks are a rounded scalar at a fixed input.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import scipy.special as sp

from harness import bench, check, require, bench_once, seed

require(["scipy.special:fresnel", "scipy.special:erfc", "scipy.special:erfi",
         "scipy.special:expi", "scipy.special:sici", "scipy.special:beta",
         "scipy.special:iv", "scipy.special:kv",
         "scipy.special:eval_legendre", "scipy.special:legendre_q"])

seed()
n = 1000000
v = np.random.random(n) + 0.5          # [0.5, 1.5]: safe for every function
w = v + 0.5
vc = np.random.random(n) - 0.5         # [-0.5, 0.5]: on the cut for LegendreQ


def r6(x):
    return int(np.floor(float(x) * 1e6 + 0.5))


def legendre_q_np(deg, x):
    """Vectorized Q_n on the cut via the 3-term recurrence (scipy has no
    vectorized legendre_q; only the scalar lqn), so this numpy recurrence is
    what a scipy user writes.  Q0 = arctanh(x), Q1 = x*Q0 - 1."""
    q0 = np.arctanh(x)
    if deg == 0:
        return q0
    q1 = x * q0 - 1.0
    for k in range(2, deg + 1):
        q0, q1 = q1, ((2 * k - 1) * x * q1 - (k - 1) * q0) / k
    return q1


# ---- vectorized on both sides: fair overhead comparison -----------------

bench("FresnelC over 10^6", lambda: sp.fresnel(v)[1])
check("FresnelC over 10^6", r6(sp.fresnel(0.5)[1]))

bench("FresnelS over 10^6", lambda: sp.fresnel(v)[0])
check("FresnelS over 10^6", r6(sp.fresnel(0.5)[0]))

bench("Erfc over 10^6", lambda: sp.erfc(v))
check("Erfc over 10^6", r6(sp.erfc(0.5)))

bench("Erfi over 10^6", lambda: sp.erfi(v))
check("Erfi over 10^6", r6(sp.erfi(0.5)))

bench("ExpIntegralEi over 10^6", lambda: sp.expi(v))
check("ExpIntegralEi over 10^6", r6(sp.expi(1.5)))

bench("SinIntegral over 10^6", lambda: sp.sici(v)[0])
check("SinIntegral over 10^6", r6(sp.sici(1.5)[0]))

bench("CosIntegral over 10^6", lambda: sp.sici(v)[1])
check("CosIntegral over 10^6", r6(sp.sici(1.5)[1]))

bench("Beta[.,.] over 10^6", lambda: sp.beta(v, w))
check("Beta[.,.] over 10^6", r6(sp.beta(1.5, 2.0)))

# ---- Legendre functions of both kinds over the cut ----------------------
# P: scipy.special.eval_legendre is a compiled vectorized ufunc (legendre_p
# exists too but returns a stacked shape).  Q: scipy has NO vectorized second
# kind, so the numpy recurrence above stands in.

bench("LegendreP[10, .] over 10^6", lambda: sp.eval_legendre(10, vc))
check("LegendreP[10, .] over 10^6", r6(sp.eval_legendre(10, 0.5)))

bench("LegendreQ[10, .] over 10^6", lambda: legendre_q_np(10, vc))
check("LegendreQ[10, .] over 10^6", r6(legendre_q_np(10, 0.5)))

# ---- scipy vectorized; Mathilda scalar-threaded (missing kernel) --------

bench_once("BesselI[0, .] over 10^6", lambda: sp.iv(0, v))
check("BesselI[0, .] over 10^6", r6(sp.iv(0, 1.5)))

bench_once("BesselK[0, .] over 10^6", lambda: sp.kv(0, v))
check("BesselK[0, .] over 10^6", r6(sp.kv(0, 1.5)))
