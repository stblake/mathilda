#!/usr/bin/env python3
"""Experiment 26 -- Structural operations: block moves and gathers (numpy column).

Same seven kernels as ``structural_ops.m``, same order and sizes.

ROADMAP ITEM 9.  Pure memory motion, no arithmetic, so this is a comparison
against memory bandwidth -- the honest ceiling.

INDEXING.  The .m file's index array is 1-based; here it is used as ``x[idx - 1]``
so both files hold the same numbers, which matters when a check disagrees and has
to be traced.

VIEWS.  ``[::-1]`` and slicing are views in numpy, so ``Reverse`` and ``Take``
call ``.copy()`` -- the Wolfram operations materialise, and timing a view against
a copy would report a meaningless win.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from harness import bench, check, require, seed
from data import lcg_ints

require(["numpy:take", "numpy:roll", "numpy:concatenate", "numpy:sort",
         "numpy:unique"])

seed()
n = 4000000
v = np.random.random(n)
idx = np.random.randint(1, n + 1, size=n)          # 1-based, as in the .m file

bench("Reverse 4x10^6", lambda: v[::-1].copy())
check("Reverse 4x10^6", int(np.arange(1, 11)[::-1].sum()))

# Wolfram's RotateLeft[v, k] == numpy's roll(v, -k).
bench("RotateLeft 4x10^6 by 1000", lambda: np.roll(v, -1000))
check("RotateLeft 4x10^6 by 1000", int(np.roll(np.arange(1, 11), -3).sum()))

bench("gather v[[idx]], 4x10^6", lambda: v[idx - 1])
check("gather v[[idx]], 4x10^6",
      int(np.arange(1, 11)[np.array([2, 4, 6]) - 1].sum()))

bench("Join two 2x10^6",
      lambda: np.concatenate([v[:2000000], v[-2000000:]]))
check("Join two 2x10^6",
      int(np.concatenate([np.arange(1, 6), np.arange(1, 6)]).sum()))

bench("Take first half of 4x10^6", lambda: v[:2000000].copy())
check("Take first half of 4x10^6", int(np.arange(1, 11)[:5].sum()))

bench("Sort 4x10^6", lambda: np.sort(v), reps=1)
check("Sort 4x10^6", int(sum(sorted(lcg_ints(2000, 50000))[:10])))

bench("Union of 4x10^6 integers", lambda: np.unique(idx), reps=1)
check("Union of 4x10^6 integers", len(np.unique(np.array([3, 1, 2, 3, 1]))))
