"""Experiment 196 -- SIAM 100-digit holdout A, family hc1 (Python column).

The ten integrals of oscillatory_tails.m, as written over the stated interval:
scipy.integrate.quad (limit 1000) at machine precision, mpmath.quad at 60
digits.  hc1.7 overflows math.exp for x < 1/709 (a raised OverflowError
scores as failed)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath as mp
from scipy.integrate import quad

require(["scipy.integrate:quad", "mpmath:quad"])

M = {1: (lambda x: math.sin(1/x), 0, 1),
     2: (lambda x: math.sin(1/x)/x, 0, 1),
     3: (lambda x: math.cos(x**-2)/x, 0, 1),
     4: (lambda x: math.cos(1/x)/math.sqrt(x), 0, 1),
     5: (lambda x: math.sin(math.log(x)/x)/x, 0, 1),
     6: (lambda x: math.cos(math.log(x)/x), 0, 1),
     7: (lambda x: math.cos(math.exp(1/x)), 0, 1),
     8: (lambda x: math.sin(1/(x*(1 - x))), 0, 1),
     9: (lambda x: math.sin(1/(2 - x))/(2 - x), 0, 2),
     10: (lambda x: math.cos(math.log(x)/x)/x**1.5, 0, 1)}
P = {1: lambda x: mp.sin(1/x),
     2: lambda x: mp.sin(1/x)/x,
     3: lambda x: mp.cos(x**-2)/x,
     4: lambda x: mp.cos(1/x)/mp.sqrt(x),
     5: lambda x: mp.sin(mp.log(x)/x)/x,
     6: lambda x: mp.cos(mp.log(x)/x),
     7: lambda x: mp.cos(mp.exp(1/x)),
     8: lambda x: mp.sin(1/(x*(1 - x))),
     9: lambda x: mp.sin(1/(2 - x))/(2 - x),
     10: lambda x: mp.cos(mp.log(x)/x)/x**1.5}

for k in range(1, 11):
    f, a, b = M[k]
    vrun("hc1.%d machine" % k, lambda: quad(f, a, b, limit=1000)[0], 10, cap=60)
    vrun("hc1.%d wp50" % k, lambda: mp.quad(P[k], [a, b]), 48, cap=90, dps=60)
