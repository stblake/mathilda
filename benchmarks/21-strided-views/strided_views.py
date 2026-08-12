#!/usr/bin/env python3
"""Experiment 21 -- Strided views: Transpose and sliding windows (numpy column).

Same six kernels as ``strided_views.m``, same order and sizes.

WHY THIS ROW SET IS EASY TO GET WRONG.  numpy's ``.T``, basic slicing and
``reshape`` all return VIEWS -- O(1), no data moved.  The Wolfram semantics
materialise a new array.  Timing ``m.T`` against ``Transpose[m]`` therefore
compares a pointer change against a memory copy and reports a meaningless
~1000x Python win.

Every row here calls ``.copy()`` (or ``np.ascontiguousarray``) wherever the
Wolfram operation produces a fresh array, so the two columns do the same amount
of work.  ``sliding_window_view`` is the one genuine exception kept as a view:
it is the *point* of roadmap item 2, so the row shows what a strided view buys
and is labelled as such rather than silently compared.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from harness import bench, check, require, seed
from data import spd_int_matrix

require(["numpy:transpose", "numpy:reshape",
         "numpy.lib.stride_tricks:sliding_window_view"])

seed()
m = np.random.random((2000, 2000))
v = np.random.random(2000000)


def r4(x):
    return int(np.floor(float(x) * 1e4 + 0.5))


# .copy() is mandatory: bare .T is a view.
bench("Transpose 2000x2000", lambda: m.T.copy())
check("Transpose 2000x2000", int(np.array([[1, 2, 3], [4, 5, 6]]).T.sum()))

bench("Transpose then Dot (fused?)", lambda: m.T @ m, reps=1)
_s6 = spd_int_matrix(6)
check("Transpose then Dot (fused?)", r4((_s6.T @ _s6).sum()))

# The strided view IS the roadmap item; materialise it so the row measures the
# same output Partition produces, and note the view cost separately in the README.
bench("Partition window 8, offset 1",
      lambda: sliding_window_view(v, 8).copy())
check("Partition window 8, offset 1",
      len(sliding_window_view(np.arange(1, 11), 8)))

bench("ArrayReshape 2x10^6 to 1000x2000",
      lambda: v.reshape(1000, 2000).copy())
check("ArrayReshape 2x10^6 to 1000x2000",
      int(np.arange(1, 7).reshape(2, 3).sum()))

bench("Take rows 1;;1000 of 2000x2000", lambda: m[:1000].copy())
check("Take rows 1;;1000 of 2000x2000", len(m[:1000]))

bench("column slice m[[All, 1]]", lambda: m[:, 0].copy())
check("column slice m[[All, 1]]", int(np.array([[1, 2], [3, 4]])[:, 0].sum()))
