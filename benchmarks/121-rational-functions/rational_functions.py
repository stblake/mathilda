#!/usr/bin/env python3
"""Experiment 121 -- Together, Cancel and Apart at stress sizes (sympy column).

The same twelve kernels as ``rational_functions.m``, in the same order, with
the same labels and the same exact integer checks.

    python3 rational_functions.py

sympy's ``together`` combines but does not cancel, so the Together cases time
``cancel``, which combines and reduces to lowest terms in one call -- the
Together contract.

COLD TIMINGS.  sympy caches results, so the harness's warm-up plus repeated
runs would time cache hits: ``cancel`` of the 50-summand sum costs 8.7 s cold
and 39 ms warm.  Every case here clears the cache before each timed run
(``cold``), with three runs when one costs under 2 s and one run otherwise.
Two cases do not finish in 120 s cold (200 summands, the depth-32 continued
fraction) and are SKIP, so the 300 s file budget holds.  The check is exact: the value at a rational point times
10^6 is a rational, rounded with no floating point, plus 10^40 times the
denominator's degree in x for Together/Cancel (an unreduced result fails it).
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import time

import sympy
from sympy import Rational as Q

from sympy.core.cache import clear_cache

from harness import check, require
from data import dense_upoly, factorable_upoly, rat_fn

require(["sympy", "sympy:cancel", "sympy:apart", "sympy:together"])

x, y, z, u, v, w, a, b, c, d, f = sympy.symbols("x y z u v w a b c d f")
pts = {x: Q(11, 10), y: Q(1, 3), z: Q(2, 7), u: Q(3, 11), v: Q(5, 13), w: Q(7, 17),
       a: Q(1, 5), b: Q(2, 9), c: Q(3, 13), d: Q(5, 19), f: Q(7, 23)}


def rnd(r):
    """Round an exact rational to the nearest integer (no float)."""
    return int(sympy.floor(r + Q(1, 2)))


def cold(label, fn):
    """Time fn with sympy's cache cleared before every run; return its result."""
    best, r = None, None
    for k in range(3):
        clear_cache()
        t0 = time.perf_counter()
        r = fn()
        dt = time.perf_counter() - t0
        best = dt if best is None else min(best, dt)
        if dt > 2.0:
            break
    print("BENCH\t%s\t%.3f" % (label, 1000.0 * best), flush=True)
    return r


def over_budget(label):
    print("SKIP\t%s\tsympy:cancel (over 120 s cold)" % label, flush=True)


def chk(t):
    _, den = sympy.fraction(sympy.together(t))
    return 10 ** 40 * sympy.degree(den, x) + rnd(10 ** 6 * t.subs(pts))


def chk_a(t, x0=None):
    return rnd(10 ** 6 * (t.subs(x, x0) if x0 is not None else t.subs(pts)))


# ---- Together -------------------------------------------------------------

tog50 = sympy.Add(*[(i ** 2 + 1) / (x - i) for i in range(1, 51)])
r = cold("Together, 50 linear summands", lambda: sympy.cancel(tog50))
check("Together, 50 linear summands", chk(r))

tog200 = sympy.Add(*[sympy.Integer(i + 1) / (x - i) for i in range(1, 201)])
over_budget("Together, 200 linear summands")

togMV = sympy.Add(*[(x + i * y) / ((x + i * y + z) * (y - i * z + 1)) for i in range(1, 13)])
r = cold("Together, 12 trivariate summands", lambda: sympy.cancel(togMV))
check("Together, 12 trivariate summands", chk(r))

nest32 = x
for _ in range(32):
    nest32 = x + 1 / (1 + nest32)
over_budget("Together, nested depth 32")

# ---- Cancel ---------------------------------------------------------------

sh200 = sympy.expand(factorable_upoly(100, x))
c400 = sympy.expand(sh200 * dense_upoly(200, x)) / sympy.expand(sh200 * (dense_upoly(199, x) + 1))
r = cold("Cancel deg-400, deg-200 gcd", lambda: sympy.cancel(c400))
check("Cancel deg-400, deg-200 gcd", chk(r))

sh500 = sympy.expand(factorable_upoly(250, x))
c1000 = sympy.expand(sh500 * dense_upoly(500, x)) / sympy.expand(sh500 * (dense_upoly(499, x) + 1))
r = cold("Cancel deg-1000, deg-500 gcd", lambda: sympy.cancel(c1000))
check("Cancel deg-1000, deg-500 gcd", chk(r))

g3 = (1 + x + y + z) ** 6
c3 = sympy.expand(g3 * (x + 2 * y - z + 3) ** 4) / sympy.expand(g3 * (x - y + 2 * z + 1) ** 4)
r = cold("Cancel dense trivariate, deg-6 gcd", lambda: sympy.cancel(c3))
check("Cancel dense trivariate, deg-6 gcd", chk(r))

g6 = 1 + x * y ** 3 + z ** 2 * u * v ** 4 + w ** 5 + x * w
c6 = (sympy.expand(g6 ** 2 * (x ** 2 + z * v + 2) * (y * u + w ** 2 + 1))
      / sympy.expand(g6 ** 2 * (x * z + u ** 3 + 3) * (v * w + y ** 2 + 1)))
r = cold("Cancel sparse 6-variable", lambda: sympy.cancel(c6))
check("Cancel sparse 6-variable", chk(r))

# ---- Apart ----------------------------------------------------------------

ap40 = rat_fn(40, x)
r = cold("Apart, 40 linear factors", lambda: sympy.apart(ap40, x))
check("Apart, 40 linear factors", chk_a(r, Q(401, 10)))

apMult = (x ** 3 + 2 * x + 5) / ((x - 1) ** 10 * (x + 2) ** 10)
r = cold("Apart, multiplicity 10 and 10", lambda: sympy.apart(apMult, x))
check("Apart, multiplicity 10 and 10", chk_a(r))

apQuad = dense_upoly(9, x) / sympy.Mul(*[x ** 2 + i for i in range(1, 7)])
r = cold("Apart, 6 irreducible quadratics", lambda: sympy.apart(apQuad, x))
check("Apart, 6 irreducible quadratics", chk_a(r))

apPar = 1 / ((x + a) * (x + b) * (x + c) * (x + d) * (x + f))
r = cold("Apart, 5 parameters", lambda: sympy.apart(apPar, x))
check("Apart, 5 parameters", chk_a(r))
