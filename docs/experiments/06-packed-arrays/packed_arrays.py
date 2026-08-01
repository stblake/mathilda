#!/usr/bin/env python3
"""Experiment 6 -- Automatic packed arrays: dense lists become machine buffers

Runs the same operations as ``packed_arrays.m``, in the same order, at the
same size.  See ``README.md``.

    python3 packed_arrays.py

NumPy has no unpacked mode to compare against -- an ndarray is always a buffer
-- so this file is the reference the two CAS are measured against rather than
a before/after in its own right.
"""

import numpy as np

import time


def bench(label, fn, reps=3):
    """One untimed warm-up, then the MINIMUM of `reps` timed runs.

    The minimum, not the mean: we are measuring the cost of the work, and every
    source of noise on a loaded machine can only add.
    """
    fn()
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    print("%-52s%s ms" % (label, round(1000.0 * min(ts), 3)))


def bench1(label, fn):
    """A single timed run, no warm-up, for kernels that take seconds."""
    t0 = time.perf_counter()
    fn()
    print("%-52s%s ms  (1 run)" % (label, round(1000.0 * (time.perf_counter() - t0), 3)))


def check(label, value):
    print("%-52scheck = %s" % (label, value))


def nest(f, x, m):
    """Nest[f, x, m] -- apply f to x, m times."""
    for _ in range(m):
        x = f(x)
    return x

# ---- total -- Total (reduction)
# The simplest reduction: one pass, no allocation.
n=10**7
x=np.random.rand(n)
y=np.random.rand(n)

# ---- accum -- Accumulate (prefix scan)
# A prefix scan -- sequential, so it cannot be threaded without changing the
# summation order.

# ---- sort -- Sort
# Order statistics: the one row here that is neither a stream nor a scan.

# ---- sin -- Sin (elementwise)
# An elementwise transcendental; libm-bound rather than memory-bound.

# ---- exp -- Exp (elementwise)

# ---- dot -- Dot (inner product)
# An inner product -- one pass, and a BLAS call once packed.

# ---- triad -- STREAM triad, a = b + 3 c
# STREAM triad: the memory-bandwidth reference for this machine.

# ---- reverse -- Reverse
# The structural family. Each of these is a pure element MOVE, so on a buffer
# it is a memcpy and on a List it is one allocation per element.

# ---- rotate -- RotateLeft

# ---- diffs -- Differences


def main():
    print("Experiment 6 -- Automatic packed arrays: dense lists become machine buffers")
    print("")
    bench("Total (reduction)", lambda: x.sum())
    bench("Accumulate (prefix scan)", lambda: np.cumsum(x))
    bench("Sort", lambda: np.sort(x))
    bench("Sin (elementwise)", lambda: np.sin(x))
    bench("Exp (elementwise)", lambda: np.exp(x))
    bench("Dot (inner product)", lambda: x @ y)
    bench("STREAM triad, a = b + 3 c", lambda: x + 3.0*y)
    bench("Reverse", lambda: x[::-1].copy())
    bench("RotateLeft", lambda: np.roll(x,-3))
    bench("Differences", lambda: np.diff(x))


if __name__ == "__main__":
    main()
