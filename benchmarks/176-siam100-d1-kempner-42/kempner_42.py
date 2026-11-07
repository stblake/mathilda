"""Experiment 176 -- SIAM 100-Digit Challenge, Appendix D problem 1 (Python).

sum 1/n over n without the substring "42": scipy has no series summation, so
both rows use mpmath.nsum (its default extrapolation: Richardson / Shanks /
Euler-Maclaurin as nsum chooses) on the indicator-weighted term -- dps 15 for
the machine row, dps 60 for wp50.  Expected to fail; see kempner_42.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require, mp
import mpmath

require(["mpmath:nsum"])


def k42(n):
    return 0 if "42" in str(int(n)) else 1 / mpmath.mpf(n)


vrun("d1.machine", lambda: float(mpmath.nsum(k42, [1, mpmath.inf])), 10, dps=15)
vrun("d1.wp50", lambda: mpmath.nsum(k42, [1, mpmath.inf]), 48, dps=60)
