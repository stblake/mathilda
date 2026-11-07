"""Experiment 207 -- SIAM 100-digit holdout A, family d2 (Python column).

The ten self-similar digit-length series of digit_length_family.m.  scipy has
no series summation, so both rows use mpmath.nsum (default acceleration) on the
direct term: dps 15 for the machine rows, dps 60 for wp50.  Per-case cap 60 s."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import mpmath
from mpmath import mpf

require(["mpmath:nsum"])


def blen(n, b=2):
    k = 0
    while n:
        n //= b
        k += 1
    return k


def f0(n):
    return mpf(n) if n <= 2 else n * f0(blen(n))


def f1(n):
    return mpf(n) if n <= 2 else mpf(n) ** 2 * f1(blen(n))


def f2(n):
    return mpf(1) if n == 1 else mpf(3) if n == 2 else n * f2(blen(n))


def g4(n):
    return mpf(n) if n <= 9 else n * blen(n, 10) * g4(blen(n, 10))


def f5(n):
    return mpf(n) if n <= 2 else n * f5(blen(n, 3)) ** 2


def f7(n):
    return mpf(1) if n == 1 else n * f7(blen(n) - 1)


def f9(n):
    return mpf(n) if n <= 2 else n * mpmath.sqrt(blen(n)) * f9(blen(n))


I = int
TERMS = {
    1: lambda n: 1 / f1(I(n)),
    2: lambda n: 1 / f2(I(n)),
    3: lambda n: (-1) ** (I(n) + 1) / f0(I(n)),
    4: lambda n: 1 / g4(I(n)),
    5: lambda n: 1 / f5(I(n)),
    6: lambda n: 1 / f0(I(n)) ** 2,
    7: lambda n: 1 / f7(I(n)),
    8: lambda n: (-1) ** (I(n) + 1) / g4(I(n)),
    9: lambda n: 1 / f9(I(n)),
    10: lambda k: 1 / f0(2 * I(k) - 1),
}

for k in range(1, 11):
    t = TERMS[k]
    vrun("hd2.%d machine" % k, lambda: float(mpmath.nsum(t, [1, mpmath.inf])), 10, cap=60, dps=15)
    vrun("hd2.%d wp50" % k, lambda: mpmath.nsum(t, [1, mpmath.inf]), 48, cap=60, dps=60)
