"""Experiment 206 -- SIAM 100-digit holdout A, family d1 (Python column).

The ten digit-restricted series of kempner_family.m.  scipy has no series
summation, so both rows use mpmath.nsum (its default extrapolation) on the
indicator-weighted term: dps 15 for the machine rows, dps 60 for wp50.
Per-case cap 60 s."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import mpmath
import numpy as np

require(["mpmath:nsum"])


def ok(s, digits):
    return 1 if s not in digits else 0


def dec(n):
    return str(int(n))


TESTS = {
    1: (lambda n: "9" not in dec(n), 1),
    2: (lambda n: "0" not in dec(n), 1),
    3: (lambda n: "00" not in dec(n), 1),
    4: (lambda n: "314" not in dec(n), 1),
    5: (lambda n: "11" not in np.base_repr(int(n), 2), 1),
    6: (lambda n: "2" not in np.base_repr(int(n), 3), 1),
    7: (lambda n: "9" not in dec(n), 2),
    8: (lambda n: all(c in "13579" for c in dec(n)), 1),
    9: (lambda n: list(dec(n)) == sorted(dec(n)), 1),
    10: (lambda n: "2026" not in dec(n), 1),
}


def term(k):
    test, s = TESTS[k]
    return lambda n: (1 / mpmath.mpf(n) ** s) if test(n) else 0


for k in range(1, 11):
    vrun("hd1.%d machine" % k, lambda: float(mpmath.nsum(term(k), [1, mpmath.inf])), 10, cap=60, dps=15)
    vrun("hd1.%d wp50" % k, lambda: mpmath.nsum(term(k), [1, mpmath.inf]), 48, cap=60, dps=60)
