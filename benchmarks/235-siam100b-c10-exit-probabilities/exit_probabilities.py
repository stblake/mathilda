"""Experiment 235 -- SIAM 100-digit HOLDOUT B for Problem 10 (Python column).

The same ten series as exit_probabilities.m.  scipy has no series summation,
so both rows use mpmath.nsum (dps 15 for the machine row, returned as a
float; dps 60 for wp50) -- as in experiment 175; the double series of
bc10.4 / bc10.10 are nsum over two index ranges, and the Bessel zeros of
bc10.8 are mpmath.besseljzero.  See exit_probabilities.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import mpmath
from mpmath import pi, nsum, inf, sech, cosh, sinh, sin, exp, sqrt, besseljzero, besselj, mpf

require(["mpmath:nsum", "mpmath:besseljzero"])

sg = lambda k: -1 if int(k) % 2 else 1
S = {1: lambda: 4 / pi * nsum(lambda k: sg(k) / (2 * k + 1) * sech(2 * (2 * k + 1) * pi), [0, inf]),
     2: lambda: 4 / pi * nsum(lambda k: sg(k) / (2 * k + 1) * cosh(3 * (2 * k + 1) * pi) / cosh(5 * (2 * k + 1) * pi), [0, inf]),
     3: lambda: 4 / pi * nsum(lambda k: sg(k) / (2 * k + 1) * sinh((2 * k + 1) * pi / 2) / sinh(2 * (2 * k + 1) * pi), [0, inf]),
     4: lambda: 16 / pi ** 2 * nsum(lambda m, n: sg(m + n) / ((2 * m + 1) * (2 * n + 1))
                                    * sech(mpf(3) / 2 * pi * sqrt((2 * m + 1) ** 2 + ((2 * n + 1) / mpf(2)) ** 2)),
                                    [0, inf], [0, inf]),
     5: lambda: 4 / pi * nsum(lambda k: sin((2 * k + 1) * pi / 3) / (2 * k + 1) * exp(-(2 * k + 1) * pi), [0, inf]),
     6: lambda: 4 / pi * nsum(lambda k: sg(k) / (2 * k + 1) * cosh(5) / cosh(5 * sqrt(1 + ((2 * k + 1) * pi) ** 2)), [0, inf]),
     7: lambda: 4 / pi * nsum(lambda k: sg(k) / (2 * k + 1) * mpf(8) ** -(2 * k + 1), [0, inf]),
     8: lambda: nsum(lambda k: 2 / (besseljzero(0, int(k)) * besselj(1, besseljzero(0, int(k))))
                     * sech(10 * besseljzero(0, int(k))), [1, inf]),
     9: lambda: 4 / pi * nsum(lambda k: sg(k) / (2 * k + 1) * sinh((2 * k + 1) * pi * mpf(99) / 100) / sinh((2 * k + 1) * pi), [0, inf]),
     10: lambda: 16 / pi ** 2 * nsum(lambda m, n: sg(m + n) / ((2 * m + 1) * (2 * n + 1))
                                     * sech(5 * pi * sqrt((2 * m + 1) ** 2 + (2 * n + 1) ** 2)), [0, inf], [0, inf])}

for k in range(1, 11):
    f = S[k]
    vrun("bc10.%d machine" % k, lambda: float(f()), 10, cap=75, dps=15)
    vrun("bc10.%d wp50" % k, f, 48, cap=75, dps=60)
