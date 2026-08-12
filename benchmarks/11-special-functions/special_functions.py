#!/usr/bin/env python3
"""Experiment 11 -- Special functions over machine arrays (scipy column).

Same kernels as ``special_functions.m``, same order and sizes.

Unlike group A this is an EXECUTION comparison: scipy.special is compiled C, so
a gap here is overhead or a missing vector kernel, not a missing algorithm.
Sizes are 10**6, above PACK_MIN_ELEMENTS, so Mathilda's buffer path is what is
being measured.

Checks are a rounded scalar from a small deterministic input -- never a sum over
the random timing data, which the three systems cannot align.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import scipy.special as sp

from harness import bench, check, require, seed

require(["scipy.special:gamma", "scipy.special:erf", "scipy.special:jv",
         "scipy.special:gammaln", "scipy.special:psi", "scipy.special:airy",
         "scipy.special:zeta"])

seed()
n = 1000000
v = np.random.random(n) + 0.5          # [0.5, 1.5]: safe for every function


def r6(x):
    return int(np.floor(x * 1e6 + 0.5))


bench("Gamma over 10^6", lambda: sp.gamma(v))
check("Gamma over 10^6", r6(sp.gamma(1.5)))

bench("Erf over 10^6", lambda: sp.erf(v))
check("Erf over 10^6", r6(sp.erf(0.5)))

bench("BesselJ[0, .] over 10^6", lambda: sp.jv(0, v))
check("BesselJ[0, .] over 10^6", r6(sp.jv(0, 0.5)))

bench("LogGamma over 10^6", lambda: sp.gammaln(v))
check("LogGamma over 10^6", r6(sp.gammaln(1.5)))

bench("PolyGamma[0, .] over 10^6", lambda: sp.psi(v))
check("PolyGamma[0, .] over 10^6", r6(sp.psi(1.5)))

bench("AiryAi over 10^6", lambda: sp.airy(v)[0])
check("AiryAi over 10^6", r6(sp.airy(0.5)[0]))

bench("Zeta over 10^6", lambda: sp.zeta(v + 1))
check("Zeta over 10^6", r6(sp.zeta(2.5)))
