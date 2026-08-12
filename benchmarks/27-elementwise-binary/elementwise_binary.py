#!/usr/bin/env python3
"""Experiment 27 -- Elementwise binary kernels (numpy column).

Same seven kernels as ``elementwise_binary.m``, same order and sizes.

ROADMAP ITEM 10.  Binary elementwise is where a system's SIMD story shows up most
plainly: the loop body is one instruction, so anything other than a vector op is
visible immediately.  ``minimum``/``maximum`` have no arithmetic at all -- a
compare and a select, which vectorises perfectly -- and numpy reaches ~48 GB/s on
them, which is the number roadmap item 10 is measured against.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from harness import bench, check, require, seed

require(["numpy:minimum", "numpy:maximum", "numpy:clip", "numpy:mod"])

seed()
n = 4000000
a = np.random.random(n)
b = np.random.random(n)


def r6(x):
    return int(np.floor(float(x) * 1e6 + 0.5))


bench("a + b over 4x10^6", lambda: a + b)
check("a + b over 4x10^6", int((np.arange(1, 6) + np.arange(1, 6)).sum()))

bench("a * b over 4x10^6", lambda: a * b)
check("a * b over 4x10^6", int((np.arange(1, 6) * np.arange(1, 6)).sum()))

bench("MapThread[Min] over 4x10^6", lambda: np.minimum(a, b))
check("MapThread[Min] over 4x10^6",
      int(np.minimum(np.array([1, 5, 3]), np.array([4, 2, 6])).sum()))

bench("MapThread[Max] over 4x10^6", lambda: np.maximum(a, b))
check("MapThread[Max] over 4x10^6",
      int(np.maximum(np.array([1, 5, 3]), np.array([4, 2, 6])).sum()))

bench("Clip to [0.25, 0.75] over 4x10^6", lambda: np.clip(a, 0.25, 0.75))
check("Clip to [0.25, 0.75] over 4x10^6",
      r6(np.clip(np.array([0.0, 0.5, 1.0]), 0.25, 0.75).sum()))

bench("a b + a over 4x10^6", lambda: a * b + a)
check("a b + a over 4x10^6",
      int((np.arange(1, 6) * np.arange(1, 6) + np.arange(1, 6)).sum()))

ia = np.random.randint(1, 1001, size=n)
ib = np.random.randint(1, 1001, size=n)
bench("integer Mod over 4x10^6", lambda: np.mod(ia, ib))
check("integer Mod over 4x10^6",
      int(np.mod(np.array([7, 8, 9]), np.array([3, 5, 4])).sum()))
