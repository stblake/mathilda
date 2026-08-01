#!/usr/bin/env python3
"""Experiment 11 -- Fourth sweep: closing the distance to NumPy.

Runs the same probe table as ``hpc_sweep_numpy_gap.m``.  See ``README.md``.

    python3 hpc_sweep_numpy_gap.py

This file is the REFERENCE the sweep was measured against: one operation per
line, at 10^6 float64.  Several rows have a NumPy spelling that returns a VIEW
rather than a new array (``v[:-1]``, ``v[::-1]``); a view is O(1) and would
flatter NumPy against two systems that both materialise, so those rows force a
copy.  The comparison is of the same work, not of the same syntax.

THE CONTROL IS THE POINT.  The last two rows -- a plain product and a log --
were already at parity before the sweep, which is what says that nothing above
them is about arithmetic.
"""

import time

import numpy as np
from scipy import signal as spsig


def bench(label, fn, reps=3):
    """One untimed warm-up, then the MINIMUM of `reps` timed runs."""
    fn()
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    print("%-52s%s ms" % (label, round(1000.0 * min(ts), 3)))


def check(label, value):
    print("%-52scheck = %s" % (label, value))

N = 10 ** 6
v = np.random.rand(N)
k5 = np.random.rand(5)
im = np.random.rand(1024, 1024)
k55 = np.random.rand(5, 5)


def ema(a, alpha=0.02):
    """The general linear recurrence has no NumPy primitive; lfilter is the
    library answer and is what a NumPy programmer would reach for."""
    return spsig.lfilter([alpha], [1.0, -(1.0 - alpha)], a)


def main():
    print("Experiment 11 -- structural operations, scans and convolution")
    print("")
    print("-- the structural family: pure element MOVES --")
    bench("v[0]      (an O(1) element read)", lambda: v[0])
    bench("v[-1]", lambda: v[-1])
    bench("v[:-1].copy()", lambda: v[:-1].copy())
    bench("v[1:].copy()", lambda: v[1:].copy())
    bench("v[250:].copy()", lambda: v[250:].copy())
    bench("v[249:-249].copy()", lambda: v[249:-249].copy())
    bench("v[::-1].copy()", lambda: v[::-1].copy())
    bench("roll(v, -3)", lambda: np.roll(v, -3))

    print("")
    print("-- scans --")
    bench("cumsum(v)", lambda: np.cumsum(v))
    bench("diff(v)", lambda: np.diff(v))
    bench("maximum.accumulate(v)", lambda: np.maximum.accumulate(v))
    bench("add.accumulate(v)", lambda: np.add.accumulate(v))
    bench("lfilter (EMA)", lambda: ema(v))

    print("")
    print("-- convolution --")
    bench("convolve(v, k5)", lambda: np.convolve(v, k5, mode="valid"))
    bench("correlate2d(im, k55)  (1024^2)",
          lambda: spsig.correlate(im, k55, mode="valid"))

    print("")
    print("-- clipping --")
    bench("clip(v, 0.2, 0.8)", lambda: np.clip(v, 0.2, 0.8))

    print("")
    print("-- THE CONTROL: the memory floor --")
    bench("v * v", lambda: v * v)
    bench("log(v)", lambda: np.log(v))


if __name__ == "__main__":
    main()
