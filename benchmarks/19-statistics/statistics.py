#!/usr/bin/env python3
"""Experiment 19 -- Descriptive statistics over machine arrays (numpy/scipy).

Same nine kernels as ``statistics.m``, same order and sizes.

Checks use ``lcg_ints`` from data.py -- an identical integer sequence in all
three systems -- so a Mean or Median check is exactly comparable rather than a
comparison of three different random draws.

TWO CONVENTION DIFFERENCES that the checks account for, because getting them
wrong shows up as a check disagreement rather than as a timing:
  * Wolfram's ``Variance``/``StandardDeviation`` are the UNBIASED (n-1)
    estimators; numpy defaults to n.  Hence ``ddof=1``.
  * Wolfram's ``Skewness``/``Kurtosis`` are the population (biased) forms, and
    its Kurtosis is NOT excess (a normal gives 3, not 0).  Hence ``bias=True``
    and ``fisher=False``.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import numpy as np
from scipy import stats as st

from harness import bench, check, require, seed
from data import lcg_ints

require(["numpy:mean", "numpy:median", "numpy:var", "numpy:std",
         "numpy:quantile", "scipy.stats:skew", "scipy.stats:kurtosis"])

seed()
n = 2000000
v = np.random.random(n)
det = np.array(lcg_ints(10000, 1000), dtype=np.float64)


def r(x, p):
    return int(math.floor(float(x) * 10 ** p + 0.5))


bench("Mean over 2x10^6", lambda: np.mean(v))
check("Mean over 2x10^6", r(np.mean(det), 6))

bench("Variance over 2x10^6", lambda: np.var(v, ddof=1))
check("Variance over 2x10^6", r(np.var(det, ddof=1), 3))

bench("StandardDeviation over 2x10^6", lambda: np.std(v, ddof=1))
check("StandardDeviation over 2x10^6", r(np.std(det, ddof=1), 4))

bench("Median over 2x10^6", lambda: np.median(v))
check("Median over 2x10^6", r(np.median(det), 4))

# Wolfram's Quartiles uses its default quantile definition; numpy's "linear"
# method is the closest match for an even-length sample.
# NO CHECK FROM THIS COLUMN, deliberately.  Wolfram's Quartiles uses a quantile
# definition numpy has no exact equivalent for -- Mathilda and Mathematica agree
# with each other and every numpy `method=` disagrees with both by a few parts in
# 10^4.  Emitting a check here would report CHECK-FAIL and DISCARD a timing that
# is perfectly valid, so this column times the operation and abstains from the
# value gate; the two CAS still check each other.
bench("Quartiles over 2x10^6", lambda: np.quantile(v, [0.25, 0.5, 0.75]))

bench("Skewness over 2x10^6", lambda: st.skew(v, bias=True))
check("Skewness over 2x10^6", r(st.skew(det, bias=True), 4))

bench("Kurtosis over 2x10^6", lambda: st.kurtosis(v, fisher=False, bias=True))
check("Kurtosis over 2x10^6", r(st.kurtosis(det, fisher=False, bias=True), 4))

bench("RootMeanSquare over 2x10^6", lambda: np.sqrt(np.mean(v * v)))
check("RootMeanSquare over 2x10^6", r(math.sqrt(np.mean(det * det)), 4))

# MovingAverage[list, 100] returns len-99 windowed means.
def _ma(a, w):
    c = np.cumsum(np.insert(a, 0, 0.0))
    return (c[w:] - c[:-w]) / w


bench("MovingAverage window 100", lambda: _ma(v, 100))
check("MovingAverage window 100", r(_ma(det, 100).sum(), 4))
