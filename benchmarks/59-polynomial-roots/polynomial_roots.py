#!/usr/bin/env python3
"""Experiment 59 -- High-degree polynomial roots (numpy column).

Same kernels as ``polynomial_roots.m``, same order and labels.

Baseline is compiled (numpy.roots -> companion-matrix eigenvalues via LAPACK).
Both sides compute the same roots, so the order-invariant check Total[|root|]
agrees to 4 places even at degree 100.  Dense coefficients follow denseUPoly
from data.m (built here in pure Python so we do not pay sympy's import cost to
reach numpy.roots); NRoots' count is len(roots).
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from harness import bench, check, require

require(["numpy:roots"])


def dense_coeffs_desc(n):
    """Descending coefficient list of denseUPoly[n]: ((k^2+3k+1) mod 17) x^k + 1."""
    c = [((k * k + 3 * k + 1) % 17) for k in range(n + 1)]
    c[0] += 1
    return c[::-1]


def r(x, p):
    return int(np.floor(float(x) * 10 ** p + 0.5))


c20 = dense_coeffs_desc(20)
c50 = dense_coeffs_desc(50)
c100 = dense_coeffs_desc(100)
c40 = dense_coeffs_desc(40)
cRootUnity = [1] + [0] * 63 + [-1]                 # x^64 - 1
cWilk = np.poly(np.arange(1, 16))                  # (x-1)...(x-15)

bench("NSolve dense degree 20", lambda: np.roots(c20))
check("NSolve dense degree 20", r(np.abs(np.roots(c20)).sum(), 4))

bench("NSolve dense degree 50", lambda: np.roots(c50))
check("NSolve dense degree 50", r(np.abs(np.roots(c50)).sum(), 4))

bench("NSolve dense degree 100", lambda: np.roots(c100))
check("NSolve dense degree 100", r(np.abs(np.roots(c100)).sum(), 4))

bench("NSolve roots of unity degree 64", lambda: np.roots(cRootUnity))
check("NSolve roots of unity degree 64", r(np.abs(np.roots(cRootUnity)).sum(), 4))

bench("NSolve Wilkinson degree 15", lambda: np.roots(cWilk))
check("NSolve Wilkinson degree 15", r(np.abs(np.roots(cWilk)).sum(), 3))

bench("NRoots dense degree 40", lambda: np.roots(c40))
check("NRoots dense degree 40", int(len(np.roots(c40))))
