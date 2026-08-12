#!/usr/bin/env python3
"""Experiment 31 -- NIntegrate difficulty ladder (scipy.integrate column).

Same rows as ``nintegrate_ladder.m``, same order, same labels.

WHAT THE BASELINE IS.  scipy's ``quad`` is QUADPACK, i.e. adaptive
Gauss-Kronrod with the same lineage as the CAS columns' default rule, so these
rows compare like with like -- unlike, say, the graph experiment, where the
baseline is pure Python and a win means less.

WHERE THE BASELINE IS TOLD THE ANSWER.  QUADPACK needs help on two families,
and giving it that help is the honest comparison, because a user of scipy would
give it too:

  * oscillatory rows use ``limit`` raised well above the default 50 subdivisions
    -- without it QUADPACK returns a convergence warning rather than a number
    at k = 1000;
  * the kinked integrand is handed ``points=[1/3]`` so the panel boundary lands
    on the kink, which is what the argument exists for.

The peak and endpoint rows are deliberately NOT helped: locating a narrow peak
and handling an integrable endpoint singularity are what the adaptive rule and
the double-exponential transform respectively claim to do unaided, so helping
them would measure the wrong thing.

High-precision rows use mpmath, since float64 cannot represent the answer at
all -- that is the point of those two rows.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import mpmath
from scipy.integrate import quad, dblquad, tplquad

from harness import bench, check, require

require(["scipy.integrate:quad", "scipy.integrate:dblquad",
         "scipy.integrate:tplquad", "mpmath:quad"])


def r(x, p):
    """Round like the .m column's Round[10^p x] -- half away from zero."""
    v = x * (10 ** p)
    return int(math.floor(v + 0.5)) if v >= 0 else -int(math.floor(-v + 0.5))


# ---- oscillation ---------------------------------------------------------
# limit=2000: at k=1000 the default 50 subdivisions cannot resolve 500 periods
# and QUADPACK reports non-convergence instead of an answer.
def osc(k):
    return quad(lambda x: math.sin(k * x), 0.0, math.pi, limit=2000)[0]


bench("NI oscillatory k=40", lambda: osc(40))
check("NI oscillatory k=40", r(osc(40), 6))

bench("NI oscillatory k=200", lambda: osc(200))
check("NI oscillatory k=200", r(osc(200), 6))

bench("NI oscillatory k=1000", lambda: osc(1000))
check("NI oscillatory k=1000", r(osc(1000), 6))

bench("NI oscillatory k=1001 nonzero", lambda: osc(1001))
check("NI oscillatory k=1001 nonzero", r(osc(1001), 6))


# ---- peak width ----------------------------------------------------------
def peak(w):
    return quad(lambda x: math.exp(-w * (x - 0.5) ** 2), 0.0, 1.0)[0]


bench("NI narrow peak w=1e3", lambda: peak(1e3))
check("NI narrow peak w=1e3", r(peak(1e3), 8))

bench("NI narrow peak w=1e5", lambda: peak(1e5))
check("NI narrow peak w=1e5", r(peak(1e5), 8))


# ---- endpoint singularity ------------------------------------------------
bench("NI endpoint x^-1/2", lambda: quad(lambda x: x ** -0.5, 0.0, 1.0)[0])
check("NI endpoint x^-1/2", r(quad(lambda x: x ** -0.5, 0.0, 1.0)[0], 6))

bench("NI endpoint x^-9/10", lambda: quad(lambda x: x ** -0.9, 0.0, 1.0)[0])
check("NI endpoint x^-9/10", r(quad(lambda x: x ** -0.9, 0.0, 1.0)[0], 4))


# ---- interior kink -------------------------------------------------------
def kink():
    return quad(lambda x: abs(x - 1.0 / 3.0), 0.0, 1.0, points=[1.0 / 3.0])[0]


bench("NI interior kink", kink)
check("NI interior kink", r(kink(), 6))


# ---- semi-infinite range -------------------------------------------------
def semiinf():
    return quad(lambda x: math.exp(-x * x), 0.0, math.inf)[0]


bench("NI semi-infinite Gaussian", semiinf)
check("NI semi-infinite Gaussian", r(semiinf(), 6))


# ---- dimension -----------------------------------------------------------
def d2():
    return dblquad(lambda y, x: x * y, 0.0, 1.0, 0.0, 1.0)[0]


def d3():
    return tplquad(lambda z, y, x: x * y * z, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)[0]


def ridge():
    return dblquad(lambda y, x: math.exp(-100.0 * (x - y) ** 2),
                   0.0, 1.0, 0.0, 1.0)[0]


bench("NI 2-D product", d2)
check("NI 2-D product", r(d2(), 6))

bench("NI 3-D product", d3)
check("NI 3-D product", r(d3(), 6))

bench("NI 2-D ridge", ridge)
check("NI 2-D ridge", r(ridge(), 6))


# ---- precision -----------------------------------------------------------
# float64 has ~15.9 digits, so these two rows are mpmath's by necessity, not by
# preference: the .m columns are being asked for 30 and 50 correct digits.
def mp_sin30():
    with mpmath.workdps(30):
        return mpmath.quad(mpmath.sin, [0, mpmath.pi])


def mp_gauss50():
    with mpmath.workdps(50):
        return mpmath.quad(lambda t: mpmath.exp(-t * t), [0, 3])


bench("NI 30-digit Sin", mp_sin30)
check("NI 30-digit Sin", r(float(mp_sin30()), 6))

bench("NI 50-digit Gaussian", mp_gauss50)
check("NI 50-digit Gaussian", r(float(mp_gauss50()), 6))
