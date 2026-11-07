"""Experiment 236 -- SIAM 100-Digit Challenge, holdout B for Appendix D
problem 1 (Python column).  Ten digit-restricted series; see
digit_restricted_series.m.  scipy has no series summation, so both rows use
mpmath.nsum (its default extrapolation) on the indicator-weighted term --
dps 15 for the machine row, dps 60 for wp50."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require, mp
import mpmath
import numpy as np

require(["mpmath:nsum"])

ok = {
    1: lambda m: "9" not in str(m),
    2: lambda m: "0" not in str(m),
    3: lambda m: all(c in "13579" for c in str(m)),
    4: lambda m: "314" not in str(m),
    5: lambda m: str(m).count("9") == 1,
    6: lambda m: "7" not in str(m),
    7: lambda m: "11" not in np.base_repr(m, 2),
    8: lambda m: "2" not in np.base_repr(m, 3),
    9: lambda m: "11" not in np.base_repr(m, 2),
    10: lambda m: all(a != b for a, b in zip(str(m), str(m)[1:])),
}


def wt(p, n):
    if p == 6:
        return 1 / mpmath.mpf(n) ** 2
    if p == 9:
        return mpmath.mpf(2) ** -n
    return 1 / mpmath.mpf(n)


def term(p):
    return lambda n: wt(p, n) if ok[p](int(n)) else 0


for p in range(1, 11):
    vrun("bd1.%d machine" % p, lambda: float(mpmath.nsum(term(p), [1, mpmath.inf])), 10,
         cap=60, dps=15)
    vrun("bd1.%d wp50" % p, lambda: mpmath.nsum(term(p), [1, mpmath.inf]), 48, cap=60, dps=60)
