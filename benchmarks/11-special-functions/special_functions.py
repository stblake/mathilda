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
         "scipy.special:zeta", "scipy.special:ellipk", "scipy.special:ellipe",
         "scipy.special:ellipkinc", "scipy.special:ellipeinc",
         "scipy.special:elliprf", "scipy.special:elliprj"])

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

# --- Legendre elliptic integrals ------------------------------------------
# SciPy's ellipk/ellipe/ellipkinc/ellipeinc take the PARAMETER m = k^2, the same
# convention as Mathilda and Mathematica, so these are direct counterparts.
#
# m in [0, 0.95] and phi in [0, 1.5]: v reaches past 1, and m > 1 is where the
# complete integrals go complex and Mathilda's machine kernels decline, which
# would time the decline path rather than the kernel.
#
# SciPy has no complete third kind, so this column composes Carlson's forms,
#   Pi(n|m) == R_F(0, 1-m, 1) + (n/3) R_J(0, 1-m, 1, 1-n),
# verified against mpmath at <= 1 ulp. There is no comparable one-liner for the
# INCOMPLETE third kind, so that row is absent here.
em = 0.95 * np.random.random(n)
ep = 1.5 * np.random.random(n)
_z = np.zeros_like(em)
_o = np.ones_like(em)


def _ellippi(nn, m):
    return sp.elliprf(_z, 1.0 - m, _o) + (nn / 3.0) * sp.elliprj(
        _z, 1.0 - m, _o, np.full_like(m, 1.0 - nn))


bench("EllipticK over 10^6", lambda: sp.ellipk(em))
check("EllipticK over 10^6", r6(sp.ellipk(0.5)))

bench("EllipticE over 10^6", lambda: sp.ellipe(em))
check("EllipticE over 10^6", r6(sp.ellipe(0.5)))

bench("EllipticF[., m] over 10^6", lambda: sp.ellipkinc(ep, em))
check("EllipticF[., m] over 10^6", r6(sp.ellipkinc(1.0, 0.5)))

bench("EllipticE[., m] over 10^6", lambda: sp.ellipeinc(ep, em))
check("EllipticE[., m] over 10^6", r6(sp.ellipeinc(1.0, 0.5)))

bench("EllipticPi[n, m] over 10^6", lambda: _ellippi(0.5, em))
check("EllipticPi[n, m] over 10^6",
      r6(float(sp.elliprf(0.0, 0.5, 1.0) + (1.0 / 3.0 / 3.0) * sp.elliprj(0.0, 0.5, 1.0, 1.0 - 1.0 / 3.0))))
