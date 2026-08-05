#!/usr/bin/env python3
"""Experiment 23 -- The per-operation constant (numpy column).

Same six kernels as ``per_operation_constant.m``, same order and sizes.

METHOD.  Hold the total element count fixed and vary the OPERATION count.  A
system with zero per-op overhead takes the same time for 1 op on 10**6 elements
as for 1000 ops on 10**3.  The gap between those rows IS the constant; divide by
the op count to read it in microseconds.

numpy's own per-call overhead is ~1-2 us (Python dispatch plus ufunc setup), so
this experiment measures Mathilda's constant against a known, documented one
rather than against zero.  The last row -- a pure scalar loop -- is CPython, not
a library call, and must be read as "an interpreted scalar loop".
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from harness import bench, check, require, seed

require(["numpy:sum", "numpy:add"])

seed()
big = np.random.random(1000000)
small = np.random.random(1000)
tiny = np.random.random(10)
hundred = np.random.random(100)


def r6(x):
    return int(np.floor(float(x) * 1e6 + 0.5))


bench("1 op on 10^6 elements", lambda: big.sum())
check("1 op on 10^6 elements", r6(np.arange(1, 11).sum()))

bench("1000 ops on 10^3 elements",
      lambda: [small.sum() for _ in range(1000)])
check("1000 ops on 10^3 elements", r6(np.arange(1, 11).sum()))

bench("100000 ops on 10 elements",
      lambda: [tiny.sum() for _ in range(100000)])
check("100000 ops on 10 elements", r6(np.arange(1, 11).sum()))

bench("1 elementwise op on 10^6", lambda: big + big)
check("1 elementwise op on 10^6",
      r6((np.arange(1, 11) + np.arange(1, 11)).sum()))

bench("10000 elementwise ops on 10^2",
      lambda: [hundred + hundred for _ in range(10000)])
check("10000 elementwise ops on 10^2",
      r6((np.arange(1, 11) + np.arange(1, 11)).sum()))


def _scalar_loop():
    s = 0.0
    for _ in range(100000):
        s = s + 1.0
    return s


bench("100000 scalar additions", _scalar_loop)
check("100000 scalar additions", 1)
