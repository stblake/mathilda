#!/usr/bin/env python3
"""Experiment 8 -- A machine-precision random number generator.

Runs the same draws as ``random_number_generation.m``.  See ``README.md``.

    python3 random_number_generation.py

NumPy's default generator is PCG64, a different algorithm from either CAS's,
so the comparison is of THROUGHPUT only -- these are not the same bits.  The
Monte Carlo row is the meaningful one: it consumes the draws, and its answer
must converge to pi in every system.
"""

import time

import numpy as np


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


def mcpi(m):
    u = np.random.rand(m)
    v = np.random.rand(m)
    return 4.0 * float((u * u + v * v <= 1.0).sum()) / m


def main():
    print("Experiment 8 -- random number generation")
    print("")
    print("-- bulk draws, 10^7 --")
    bench("rand(10^7)", lambda: np.random.rand(10 ** 7))
    bench("randint(0, 101, 10^7)", lambda: np.random.randint(0, 101, 10 ** 7))
    bench("randint(0, 256, 10^7)", lambda: np.random.randint(0, 256, 10 ** 7))
    bench("rand(1000, 1000)", lambda: np.random.rand(1000, 1000))

    print("")
    print("-- a Monte Carlo that consumes them --")
    bench("Monte Carlo pi, 10^7 samples", lambda: mcpi(10 ** 7))
    print("%-52svalue = %s   (converges to pi, not exact)"
          % ("Monte Carlo pi, 10^7 samples", mcpi(10 ** 7)))

    print("")
    print("-- the draw is only useful if the distribution is right --")
    sq = np.random.rand(10 ** 6)
    print("mean     (want 0.5)       =", round(float(sq.mean()), 4))
    print("variance (want 1/12)      =", round(float(sq.var()), 4))
    si = np.random.randint(0, 10, 10 ** 6)
    print("integer mean (want 4.5)   =", round(float(si.mean()), 4))


if __name__ == "__main__":
    main()
