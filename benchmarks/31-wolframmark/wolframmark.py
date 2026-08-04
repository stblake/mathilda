#!/usr/bin/env python3
"""Experiment 31 -- WolframMark, Python column.

The same fifteen tests as ``wolframmark.m``, at the same sizes and repetition
counts, in the same order, using numpy / scipy / sympy / mpmath.

    python3 wolframmark.py

WHAT THIS SUITE ACTUALLY IS.  WolframMark is Wolfram's HARDWARE benchmark: it
holds Mathematica constant and varies the machine, scoring against a reference
system.  It answers "how fast is this laptop at running Mathematica", not "how
good is this CAS".  So it is not a coverage measure, its sizes are tuned to
Mathematica's performance profile rather than to the workloads' intrinsic cost
(the 1.2M-point x11 DFT is a ~1 s test in Mathematica and would have taken ~40
hours against Mathilda's pre-FFTW O(n^2) fallback), and its aggregate SCORE is
meaningless across implementations -- so no score is computed here, just 15 more
cases.

What it IS good for: the workload selection is a third party's and predates this
project, so unlike groups A-C it cannot be accused of being picked to flatter us.

THREE PLACES THE PYTHON COLUMN MUST BE WRITTEN CAREFULLY, or it wins by not
doing the work:

1. **Matrix Transpose.**  ``m.T`` in numpy is a *view* -- O(1), no data moved.
   The Wolfram test materialises a transposed matrix, so timing ``m.T`` would
   compare a memory copy against a pointer change and report a ~1000x Python
   win that means nothing.  This uses ``m.T.copy()``.

2. **Gamma of large integers.**  ``scipy.special.gamma(85000)`` overflows to
   ``inf`` in float64 and returns instantly.  Mathematica's ``Gamma[n]`` for
   integer n is the *exact* integer (n-1)!, hundreds of thousands of digits
   long.  The honest equivalent is ``math.factorial(n - 1)``.

3. **Digits of Pi.**  numpy has no arbitrary precision; ``mpmath`` is the
   equivalent, and it is a pure-Python bignum library, so this row reads as
   "mpmath", not as "what C can do".

Where numpy would not be the honest baseline the docstring says so, rather than
letting the table imply a library result.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import numpy as np

from harness import bench, bench_once, check, require, seed
from data import hilbert_like, lcg_ints, spd_int_matrix

require(["numpy", "scipy:optimize", "scipy.optimize:curve_fit", "sympy",
         "mpmath", "scipy.integrate:dblquad", "numpy.fft:fft",
         "numpy.linalg:eigvals", "numpy.linalg:svd", "numpy.linalg:solve"])

# ---- 0. The compatibility bug, as its own row ---------------------------

# Mathilda leaves `RandomReal[{}, dims]` unevaluated where Mathematica accepts
# it.  Python has no such form; the check simply reports True, so the row's
# disagreement with Mathilda is what surfaces the bug.
bench("RandomReal empty-range form", lambda: np.random.random(200000))
check("RandomReal empty-range form", True)

# ---- 1. Data Fitting ----------------------------------------------------


def _fit():
    from scipy.optimize import curve_fit
    ax = np.arange(0.2, 10.0 + 1e-12, 0.22)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    x, y, z = X.ravel(), Y.ravel(), Z.ravel()
    d = np.log(120 * x) - np.abs(np.cos(z / 300) / (140 * y))

    def model(xyz, a, b, c):
        xx, yy, zz = xyz
        return np.log(a * xx) - np.abs(np.cos(b * zz) / (c * yy))

    curve_fit(model, (x, y, z), d, p0=[100.0, 1.0 / 300.0, 100.0], maxfev=20000)


bench_once("Data Fitting", _fit)

# ---- 2. Digits of Pi ----------------------------------------------------


def _pi():
    import mpmath
    mpmath.mp.dps = 1000000
    return +mpmath.pi


bench_once("Digits of Pi", _pi)


def _pi_check():
    import mpmath
    mpmath.mp.dps = 60
    return int(mpmath.floor(mpmath.mpf(10) ** 20 * (mpmath.pi - 3)
                            + mpmath.mpf("0.5")))


check("Digits of Pi", _pi_check())

# ---- 3. Discrete Fourier Transform -------------------------------------

# Wolfram's Fourier[] normalises by 1/Sqrt[n]; numpy's fft does not.  Only the
# CHECK is affected (the timing is the same transform), so the check scales.
seed()
_dft_data = np.random.random(1200000)
bench_once("Discrete Fourier Transform",
           lambda: [np.fft.fft(_dft_data) for _ in range(11)])
_small = np.array([1., 2., 3., 4., 5., 6., 7., 8.])
check("Discrete Fourier Transform",
      round(1e6 * (np.fft.fft(_small)[0].real / math.sqrt(8))))

# ---- 4. Eigenvalues of a Matrix ----------------------------------------

# a**-1 is the elementwise reciprocal in the Wolfram source (Power is Listable),
# so it is elementwise here too -- copying the official test exactly.
seed()
_a = np.random.random((420, 420))
_b = np.diag(np.random.random(420))
_m = _a @ _b @ (_a ** -1)
bench_once("Eigenvalues of a Matrix",
           lambda: [np.linalg.eigvals(_m) for _ in range(6)])
check("Eigenvalues of a Matrix",
      round(1e6 * float(np.sort(np.linalg.eigvals(hilbert_like(6)).real).sum())))

# ---- 5. Elementary Functions -------------------------------------------

seed()
_m1 = np.random.random(2200000)
_m2 = np.random.random(2200000)


def _elem():
    for _ in range(30):
        np.exp(_m1)
        np.sin(_m1)
        np.arctan2(_m1, _m2)


bench_once("Elementary Functions", _elem)
# Wolfram ArcTan[x, y] == atan2(y, x): reversed argument order.
check("Elementary Functions",
      round(1e6 * (math.exp(0.5) + math.sin(0.5) + math.atan2(1.5, 0.5))))

# ---- 6. Gamma Function -------------------------------------------------

# EXACT integers, matching Gamma[n] == (n-1)!.  scipy.special.gamma would
# overflow to inf instantly and measure nothing.
seed()
_gints = np.random.randint(80000, 90001, size=55)
bench_once("Gamma Function",
           lambda: [math.factorial(int(n) - 1) for n in _gints])
check("Gamma Function", math.factorial(20))

# ---- 7. Large Integer Multiplication -----------------------------------

# Python's int is arbitrary precision natively (CPython uses schoolbook /
# Karatsuba, not GMP -- so this row is "CPython bignum", not "GMP").
_bigint = 10 ** 1100000 + 7
bench_once("Large Integer Multiplication",
           lambda: [_bigint * (_bigint + 1) for _ in range(20)])
check("Large Integer Multiplication", ((10 ** 50 + 1) * (10 ** 50 + 2)) % 10 ** 12)

# ---- 8. Matrix Arithmetic ----------------------------------------------

seed()
_ma = np.random.random((840, 840))
bench_once("Matrix Arithmetic",
           lambda: [(1.0 + 0.5 * _ma) ** 127 for _ in range(50)])
check("Matrix Arithmetic", round(1e6 * ((1.0 + 0.5 * 0.25) ** 127)))

# ---- 9. Matrix Multiplication ------------------------------------------

seed()
_mm1 = np.random.random((1050, 1050))
_mm2 = np.random.random((1050, 1050))
bench_once("Matrix Multiplication", lambda: [_mm1 @ _mm2 for _ in range(12)])
_s8 = spd_int_matrix(8)
check("Matrix Multiplication", round(1e6 * float((_s8 @ _s8).sum())))

# ---- 10. Matrix Transpose ----------------------------------------------

# .copy() is mandatory: bare .T is a view and would time a pointer change.
seed()
_mt = np.random.random((2070, 2070))
bench_once("Matrix Transpose", lambda: [_mt.T.copy() for _ in range(40)])
check("Matrix Transpose", int(np.array([[1, 2, 3], [4, 5, 6]]).T.sum()))

# ---- 11. Numerical Integration -----------------------------------------


def _nint():
    from scipy.integrate import dblquad
    lo, hi = -2.6 * math.pi, 2.6 * math.pi
    dblquad(lambda y, x: math.sin(x * x + y * y), lo, hi, lambda _: lo,
            lambda _: hi)


bench_once("Numerical Integration", _nint)


def _nint_check():
    from scipy.integrate import quad
    return round(1e4 * quad(math.sin, 0.0, math.pi)[0])


check("Numerical Integration", _nint_check())

# ---- 12. Polynomial Expansion ------------------------------------------


def _expand():
    import sympy
    x = sympy.Symbol("x")
    p = sympy.Integer(1)
    for c in range(1, 351):
        p = p * (c + x) ** 3
    sympy.expand(p)


bench_once("Polynomial Expansion", _expand)


def _expand_check():
    import sympy
    x = sympy.Symbol("x")
    p = sympy.Integer(1)
    for c in range(1, 13):
        p = p * (c + x) ** 3
    # Length[] of an expanded Plus is its term count.
    return len(sympy.Add.make_args(sympy.expand(p)))


check("Polynomial Expansion", _expand_check())

# ---- 13. Random Number Sort --------------------------------------------

seed()
_sortme = np.random.randint(1, 50001, size=520000)
bench_once("Random Number Sort",
           lambda: [np.sort(_sortme) for _ in range(15)])
check("Random Number Sort", int(sum(sorted(lcg_ints(2000, 50000))[:10])))

# ---- 14. Singular Value Decomposition ----------------------------------

seed()
_svdm = np.random.random((860, 860))
bench_once("Singular Value Decomposition",
           lambda: [np.linalg.svd(_svdm) for _ in range(2)])
# Wolfram's SVD returns {u, sigma, v} with sigma a diagonal MATRIX; part 2 is
# that matrix, whose total is the sum of the singular values.
check("Singular Value Decomposition",
      round(1e4 * float(np.linalg.svd(spd_int_matrix(6), compute_uv=False).sum())))

# ---- 15. Solving a Linear System ---------------------------------------

seed()
_lsm = np.random.random((1150, 1150))
_lsv = np.random.random(1150)
bench_once("Solving a Linear System",
           lambda: [np.linalg.solve(_lsm, _lsv) for _ in range(16)])
check("Solving a Linear System",
      round(1e6 * float(np.linalg.solve(spd_int_matrix(8),
                                       np.arange(1.0, 9.0)).sum())))
