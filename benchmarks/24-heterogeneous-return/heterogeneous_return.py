#!/usr/bin/env python3
"""Experiment 24 -- The heterogeneous return (numpy column).

Same five kernels as ``heterogeneous_return.m``, same order and sizes.

ROADMAP ITEM 1.  A function returning ``(reals, ints, mask)`` costs Mathilda 96x
at the return and 120x on the caller's next operation, because a mixed-dtype list
of packed rows cannot be absorbed and every element is materialised as a boxed
Expr.

numpy has no such problem: a tuple of arrays is a tuple of independent buffers,
each keeping its own dtype.  That is exactly the container roadmap item 1 asks
for, which makes numpy the right reference for what the item is worth -- the
ratio on these rows IS the value of implementing it.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from harness import bench, check, require, seed

require(["numpy:round", "numpy:heaviside"])

seed()
n = 500000
v = np.random.random(n)


def r4(x):
    return int(np.floor(float(x) * 1e4 + 0.5))


def homog(x):
    return (x + 1.0, x * 2.0)


def mixed2(x):
    return (x + 1.0, np.rint(x * 10).astype(np.int64))


def mixed3(x):
    return (x + 1.0, np.rint(x * 10).astype(np.int64),
            np.heaviside(x - 0.5, 1.0).astype(np.int64))


def ragged(x):
    return (x, x[:1000].copy(), np.rint(x[:100] * 10).astype(np.int64))


bench("return {real, real}, then Total",
      lambda: sum(a.sum() for a in homog(v)))
check("return {real, real}, then Total",
      r4(homog(np.arange(1.0, 11.0))[0].sum()))

bench("return {real, int}, then Total",
      lambda: sum(a.sum() for a in mixed2(v)))
check("return {real, int}, then Total",
      r4(mixed2(np.arange(1.0, 11.0) / 10.0)[1].sum()))

bench("return {real, int, mask}, then Total",
      lambda: sum(a.sum() for a in mixed3(v)))
check("return {real, int, mask}, then Total",
      r4(mixed3(np.arange(1.0, 11.0) / 10.0)[2].sum()))

bench("return {real, int, mask}, discarded", lambda: mixed3(v))
check("return {real, int, mask}, discarded",
      len(mixed3(np.arange(1.0, 11.0) / 10.0)))

bench("return ragged {n, 1000, 100}, then Total",
      lambda: sum(a.sum() for a in ragged(v)))
check("return ragged {n, 1000, 100}, then Total",
      len(ragged(np.random.random(2000))))
