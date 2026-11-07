"""Experiment 177 -- SIAM 100-Digit Challenge, Appendix D Problem 2 (Python column).

scipy has no infinite-series summation, so both rows use mpmath.nsum, the
idiomatic Python numerical series summer (default acceleration: Richardson /
Shanks), on the direct term: at mp.dps = 15 for the machine row and at
dps 60 for the 50-digit row.  See digit_length_series.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import mpmath

require(["mpmath:nsum"])


def f(n):
    n = int(n)
    return n if n <= 2 else n * f(n.bit_length())


term = lambda n: mpmath.mpf(1) / f(n)
vrun("d2.machine", lambda: float(mpmath.nsum(term, [1, mpmath.inf])), 10, 3, dps=15)
vrun("d2.wp50", lambda: mpmath.nsum(term, [1, mpmath.inf]), 48, dps=60)
