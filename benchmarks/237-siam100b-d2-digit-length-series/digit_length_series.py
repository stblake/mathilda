"""Experiment 237 -- SIAM 100-Digit Challenge, holdout B for Appendix D
problem 2 (Python column).  Ten digit-length / logarithmic series; see
digit_length_series.m.  scipy has no series summation, so both rows use
mpmath.nsum (default acceleration) on the direct term: dps 15 for the
machine row, dps 60 for wp50."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import mpmath

require(["mpmath:nsum"])


def ilen(m, b=10):
    c = 0
    while m:
        m //= b
        c += 1
    return c


def g5(m):
    return 1 if m == 1 else m * g5(m.bit_length() - 1)


def f6(m):
    return 1 if m == 1 else m * f6(ilen(m)) ** 2


def f10(m):
    return m if m <= 2 else 2 * m * f10(ilen(m, 3))


one = lambda: mpmath.mpf(1)
tm = {
    1: lambda m: one() / (m * ilen(m) ** 2),
    2: lambda m: one() / (m * m.bit_length() ** 2),
    5: lambda m: one() / g5(m),
    6: lambda m: one() / f6(m),
    7: lambda m: one() / (m * ilen(m) * ilen(ilen(m)) ** 2),
    8: lambda m: one() / (m * mpmath.factorial(ilen(m))),
    9: lambda m: (-1) ** m.bit_length() * one() / (m * m.bit_length()),
    10: lambda m: one() / f10(m),
}
term = lambda p: (lambda n: tm[p](int(n)))
f3 = lambda n: 1 / (n * mpmath.log(n) ** 2)
f4 = lambda n: 1 / (n * mpmath.log(n) * mpmath.log(mpmath.log(n)) ** 2)


def case(p, fn, a):
    vrun("bd2.%d machine" % p, lambda: float(mpmath.nsum(fn, [a, mpmath.inf])), 10, cap=60, dps=15)
    vrun("bd2.%d wp50" % p, lambda: mpmath.nsum(fn, [a, mpmath.inf]), 48, cap=60, dps=60)


case(1, term(1), 1)
case(2, term(2), 1)
case(3, f3, 2)
case(4, f4, 3)
for p in (5, 6, 7, 8, 9, 10):
    case(p, term(p), 1)
