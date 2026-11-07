"""Experiment 205 -- SIAM 100-Digit Challenge, holdout A for Problem 10 (Python column).

The same ten separated series as brownian_hitting.m.  scipy has no series
summation, so (as in experiment 175) every row is mpmath.nsum: dps 15 and
float() for the machine rows, dps 60 for wp50; the double series of hc10.8 and
hc10.9 are a 2-D mpmath.nsum.  See brownian_hitting.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require, mp
import mpmath
from mpmath import nsum, inf, exp, cos, sin, sech, sqrt

require(["mpmath:nsum"])

pi = lambda: mpmath.pi
sg = lambda k: (-1) ** int(k)


def rc(c, x, a):
    return (exp(-c * (a - x)) + exp(-c * (a + x))) / (1 + exp(-2 * c * a))


P = {
    1: lambda: 4 / pi() * nsum(lambda k: sg(k) / (2 * k + 1) * sech(2 * (2 * k + 1) * pi()), [0, inf]),
    2: lambda: 4 / pi() * nsum(lambda k: sg(k) / (2 * k + 1) * rc((2 * k + 1) * pi(), 3, 5), [0, inf]),
    3: lambda: 4 / pi() * nsum(lambda k: sg(k) / (2 * k + 1) * cos((2 * k + 1) * pi() / 4)
                               * rc((2 * k + 1) * pi(), 4, 5), [0, inf]),
    4: lambda: 4 / pi() * nsum(lambda k: sg(k) / (2 * k + 1) * sech((2 * k + 1) * pi() / 4), [0, inf]),
    5: lambda: 4 / pi() * nsum(lambda j: sin((2 * j + 1) * pi() / 2) / (2 * j + 1)
                               * (exp(-(2 * j + 1) * pi() / 4) - exp(-7 * (2 * j + 1) * pi() / 4))
                               / (1 - exp(-2 * (2 * j + 1) * pi())), [0, inf]),
    6: lambda: 4 / pi() * nsum(lambda k: sg(k) * exp(-(2 * k + 1) * pi()) / (2 * k + 1), [0, inf]),
    7: lambda: mpmath.mpf(1) / 4 - 8 / pi() ** 3 * nsum(lambda k: sg(k) / (2 * k + 1) ** 3
                                                     * sech((2 * k + 1) * pi()), [0, inf]),
    8: lambda: 16 / pi() ** 2 * nsum(lambda j, k: sg(j + k) / ((2 * j + 1) * (2 * k + 1))
                                     * sech(1.5 * pi() * sqrt((2 * j + 1) ** 2 + (2 * k + 1) ** 2)),
                                     [0, inf], [0, inf]),
    9: lambda: 16 / pi() ** 2 * nsum(lambda j, k: sg(j + k) / ((2 * j + 1) * (2 * k + 1))
                                     * (lambda r: (exp(-r / 4) - exp(-7 * r / 4)) / (1 - exp(-2 * r)))(
                                         pi() * sqrt((2 * j + 1) ** 2 + (2 * k + 1) ** 2)),
                                     [0, inf], [0, inf]),
    10: lambda: 4 / pi() * nsum(lambda k: sg(k) / (2 * k + 1)
                                * rc((2 * k + 1) * pi(), mpmath.mpf(19) / 10, 2), [0, inf]),
}

for i in range(1, 11):
    f = P[i]
    vrun("hc10.%d machine" % i, lambda f=f: float(f()), 10, 3, dps=15)
    vrun("hc10.%d wp50" % i, f, 48, dps=60)
