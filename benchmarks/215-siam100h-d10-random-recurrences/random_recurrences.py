"""Experiment 215 -- SIAM 100-digit HOLDOUT A, family D10 (Python column).

The Monte-Carlo statement of random_recurrences.m (no library computes a
random-product Lyapunov exponent): machine: n = 10^6 choices from
numpy.random.default_rng(1), float ratio recursion r <- c + s/r,
exp(mean(log r)); wp50: exact Python integers / Fractions x_n for n = 10^4,
then |x_n|^(1/n) in mpmath at 60 digits."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
from fractions import Fraction
import math
import numpy as np
import mpmath

require(["numpy.random:default_rng", "mpmath:root"])
nM, n5 = 10 ** 6, 10 ** 4

def draws(ch, n, p=None):
    return np.random.default_rng(1).choice(ch, n, p=p).tolist()
def lyap_m(f, ch, p=None):
    r, acc = 1.0, 0.0                  # log r_0 = 0
    for e in draws(ch, nM - 1, p):
        r = f(r, e)
        acc += math.log(r)
    return math.exp(acc / nM)
def lyap_50(g, ch, p=None):
    a, b = 1, 1
    for e in draws(ch, n5 - 1, p):
        a, b = b, g(a, b, e)
    return mpmath.root(abs(mpmath.mpf(b.numerator) / b.denominator) if isinstance(b, Fraction)
                       else abs(mpmath.mpf(b)), n5)
def first_m():
    return math.exp(sum(math.log(c) for c in draws([1.0, 2.0], nM)) / nM)
def first_50():
    x = 1
    for c in draws([1, 2], n5):
        x *= c
    return mpmath.root(mpmath.mpf(x), n5)
h = Fraction(5, 2)
C = {  # k: (ratio map, exact step, choices, probabilities)
  2: (lambda r, e: 3.0 + e / r,   lambda a, b, e: 3 * b + e * a, [-1, 1], None),
  3: (lambda r, e: 4.0 + e / r,   lambda a, b, e: 4 * b + e * a, [-1, 1], None),
  4: (lambda r, e: 2.5 + e / r,   lambda a, b, e: h * b + e * a, [-1, 1], None),
  5: (lambda r, e: 3.0 + 2 * e / r, lambda a, b, e: 3 * b + 2 * e * a, [-1, 1], None),
  6: (lambda r, e: e + 1 / r,     lambda a, b, e: e * b + a, [1, 2], None),
  7: (lambda r, e: 3.0 + e / r,   lambda a, b, e: 3 * b + e * a, [-1, 0, 1], None),
  8: (lambda r, e: 3.0 + e / r,   lambda a, b, e: 3 * b + e * a, [-1, 1], [0.25, 0.75]),
  9: (lambda r, e: e + 1 / r,     lambda a, b, e: e * b + a, [1, 2, 3], None),
  10: (lambda r, e: 2.0 + e / r,  lambda a, b, e: 2 * b + e * a, [-1, 1], [0.25, 0.75]),
}
vrun("hd10.1 machine", first_m, 2, 3, cap=60)
vrun("hd10.1 wp50", first_50, 1, 3, cap=60, dps=60)
for k in range(2, 11):
    f, g, ch, p = C[k]
    vrun("hd10.%d machine" % k, lambda: lyap_m(f, [float(c) for c in ch] if k in (6, 9) else ch, p), 2, 3, cap=60)
    vrun("hd10.%d wp50" % k, lambda: lyap_50(g, ch, p), 1, 3, cap=60, dps=60)
